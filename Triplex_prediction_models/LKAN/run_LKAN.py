# Purpose: Execution script for training and evaluating LKAN Baseline on H-DNA dataset.
# Strategy: KANLinear-based Squeeze-and-Excitation + Inverted Residual blocks + BCELoss.

# Execution Example:
#   python run_LKAN.py --data_dir ./data --output_dir ./output --batch_size 512 --lr 5e-4 --epochs 10

import os
import math
import argparse
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
from collections import OrderedDict
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score, precision_score, recall_score, matthews_corrcoef

from KANLinear import KANLinear

# Standard nucleotide encoding map
encoding_map = {
    "A": [1, 0, 0, 0],
    "T": [0, 1, 0, 0],
    "C": [0, 0, 1, 0],
    "G": [0, 0, 0, 1]
}


def encode_sequence(sequence):
    encoded_seq = [encoding_map[base] for base in sequence if base in encoding_map]
    return torch.tensor(encoded_seq, dtype=torch.float).t()


def pad_sequences(sequences, target_length):
    padded_sequences = []
    for seq in sequences:
        seq_len = seq.size(1)
        if seq_len < target_length:
            left_padding = (target_length - seq_len) // 2
            right_padding = target_length - seq_len - left_padding
            padded_seq = F.pad(seq, (left_padding, right_padding))
        else:
            padded_seq = seq[:, :target_length]
        padded_sequences.append(padded_seq)
    return torch.stack(padded_sequences)


def collate_fn_factory(target_length):
    def collate_fn(batch):
        sequences, labels = zip(*batch)
        padded_sequences = pad_sequences(sequences, target_length)
        labels = torch.stack(labels)
        return padded_sequences, labels
    return collate_fn


class DNASequencesDataset(Dataset):
    def __init__(self, dataframe):
        self.data = dataframe[~dataframe['sequence'].str.contains('N')].reset_index(drop=True)
        self.encoded_sequences = [encode_sequence(seq) for seq in self.data['sequence']]
        self.labels = torch.tensor(self.data['label'].values, dtype=torch.long)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return self.encoded_sequences[idx], self.labels[idx]


class SELayer(nn.Module):
    def __init__(self, inp, oup, reduction=4):
        super(SELayer, self).__init__()
        self.avg_pool = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Sequential(
            KANLinear(oup, inp // reduction),
            KANLinear(inp // reduction, oup),
            nn.Sigmoid()
        )

    def forward(self, x):
        b, c, _ = x.size()
        y = x.view(b, c, -1).mean(dim=2)
        y = self.fc(y).view(b, c, 1)
        return x * y


class SeqNN(nn.Module):
    __constants__ = ('resize_factor',)

    def __init__(self, seqsize, use_single_channel=False, use_reverse_channel=False, use_multisubstate_channel=False,
                 block_sizes=[256, 256, 128, 128, 64, 64, 32, 32], ks=5, resize_factor=4, activation=nn.SiLU,
                 filter_per_group=2, se_reduction=4, final_ch=1, bn_momentum=0.1):
        super().__init__()
        self.block_sizes = block_sizes
        self.resize_factor = resize_factor
        self.se_reduction = se_reduction
        self.seqsize = seqsize
        self.final_ch = final_ch
        self.bn_momentum = bn_momentum
        seqextblocks = OrderedDict()

        in_channels_first_block = 4
        if use_single_channel: in_channels_first_block += 1
        if use_reverse_channel: in_channels_first_block += 1
        if use_multisubstate_channel: in_channels_first_block += 1

        seqextblocks['blc0'] = nn.Sequential(
            nn.Conv1d(in_channels=in_channels_first_block, out_channels=block_sizes[0], kernel_size=ks, padding='same', bias=False),
            nn.BatchNorm1d(block_sizes[0], momentum=self.bn_momentum),
            activation()
        )

        for ind, (prev_sz, sz) in enumerate(zip(block_sizes[:-1], block_sizes[1:])):
            seqextblocks[f'inv_res_blc{ind}'] = nn.Sequential(
                nn.Conv1d(prev_sz, sz * self.resize_factor, kernel_size=1, padding='same', bias=False),
                nn.BatchNorm1d(sz * self.resize_factor, momentum=self.bn_momentum),
                activation(),
                nn.Conv1d(sz * self.resize_factor, sz * self.resize_factor, kernel_size=ks,
                          groups=sz * self.resize_factor // filter_per_group, padding='same', bias=False),
                nn.BatchNorm1d(sz * self.resize_factor, momentum=self.bn_momentum),
                activation(),
                SELayer(prev_sz, sz * self.resize_factor, reduction=self.se_reduction),
                nn.Conv1d(sz * self.resize_factor, prev_sz, kernel_size=1, padding='same', bias=False),
                nn.BatchNorm1d(prev_sz, momentum=self.bn_momentum),
                activation(),
            )
            seqextblocks[f'resize_blc{ind}'] = nn.Sequential(
                nn.Conv1d(2 * prev_sz, sz, kernel_size=ks, padding='same', bias=False),
                nn.BatchNorm1d(sz, momentum=self.bn_momentum),
                activation()
            )

        self.seqextractor = nn.ModuleDict(seqextblocks)
        self.mapper = nn.Sequential(
            nn.Conv1d(block_sizes[-1], self.final_ch, kernel_size=1, padding='same'),
            activation()
        )

    def feature_extractor(self, x):
        x = self.seqextractor['blc0'](x)
        for i in range(len(self.block_sizes) - 1):
            x = torch.cat([x, self.seqextractor[f'inv_res_blc{i}'](x)], dim=1)
            x = self.seqextractor[f'resize_blc{i}'](x)
        return x

    def forward(self, x):
        f = self.feature_extractor(x)
        x = self.mapper(f)
        x = F.adaptive_avg_pool1d(x, 1).squeeze(2)
        return torch.sigmoid(x).squeeze(1)


def train(model, loader, criterion, optimizer, device):
    model.train()
    running_loss, all_preds, all_labels = 0.0, [], []

    for sequences, labels in loader:
        sequences, labels = sequences.to(device), labels.to(device)
        optimizer.zero_grad()
        probs = model(sequences)
        loss = criterion(probs, labels.float())
        loss.backward()
        optimizer.step()

        running_loss += loss.item() * sequences.size(0)
        all_preds.extend(probs.detach().cpu().numpy())
        all_labels.extend(labels.cpu().numpy())

    preds = (np.array(all_preds) >= 0.5).astype(int)
    return (running_loss / len(loader.dataset),
            accuracy_score(all_labels, preds),
            f1_score(all_labels, preds, zero_division=0),
            roc_auc_score(all_labels, all_preds),
            precision_score(all_labels, preds, zero_division=0),
            recall_score(all_labels, preds, zero_division=0),
            matthews_corrcoef(all_labels, preds))


def evaluate(model, loader, criterion, device):
    model.eval()
    eval_loss, all_preds, all_labels = 0.0, [], []

    with torch.no_grad():
        for sequences, labels in loader:
            sequences, labels = sequences.to(device), labels.to(device)
            probs = model(sequences)
            loss = criterion(probs, labels.float())
            eval_loss += loss.item() * sequences.size(0)

            all_preds.extend(probs.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())

    preds = (np.array(all_preds) >= 0.5).astype(int)
    return (eval_loss / len(loader.dataset),
            accuracy_score(all_labels, preds),
            f1_score(all_labels, preds, zero_division=0),
            roc_auc_score(all_labels, all_preds),
            precision_score(all_labels, preds, zero_division=0),
            recall_score(all_labels, preds, zero_division=0),
            matthews_corrcoef(all_labels, preds))


def parse_args():
    parser = argparse.ArgumentParser(description="LKAN Baseline Execution for DNA Prediction")
    parser.add_argument("--data_dir", type=str, default="./data", help="Path to directory containing train.csv, dev.csv, and test.csv")
    parser.add_argument("--output_dir", type=str, default="./output", help="Directory path to save predictions and checkpoints")
    parser.add_argument("--target_seq_len", type=int, default=512, help="Sequence length target for padding/truncation")
    parser.add_argument("--batch_size", type=int, default=512, help="Batch size for model training")
    parser.add_argument("--lr", type=float, default=5e-4, help="Learning rate for AdamW optimizer")
    parser.add_argument("--epochs", type=int, default=10, help="Total number of training epochs")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    os.makedirs(args.output_dir, exist_ok=True)
    best_model_path = os.path.join(args.output_dir, "best_kan_baseline.pth")

    train_dataset = DNASequencesDataset(pd.read_csv(os.path.join(args.data_dir, "train.csv")))
    dev_dataset = DNASequencesDataset(pd.read_csv(os.path.join(args.data_dir, "dev.csv")))
    test_dataset = DNASequencesDataset(pd.read_csv(os.path.join(args.data_dir, "test.csv")))

    collate_fn = collate_fn_factory(args.target_seq_len)

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, collate_fn=collate_fn, shuffle=True)
    dev_loader = DataLoader(dev_dataset, batch_size=args.batch_size, collate_fn=collate_fn, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=args.batch_size, collate_fn=collate_fn, shuffle=False)

    model = SeqNN(seqsize=args.target_seq_len, final_ch=1).to(device)
    criterion = nn.BCELoss()
    optimizer = optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.005)

    best_dev_f1 = 0.0

    for epoch in range(args.epochs):
        train_loss, train_acc, train_f1, train_auc, train_prec, train_rec, train_mcc = train(
            model, train_loader, criterion, optimizer, device
        )
        val_loss, val_acc, val_f1, val_auc, val_prec, val_rec, val_mcc = evaluate(
            model, dev_loader, criterion, device
        )

        print(f"Epoch {epoch + 1}/{args.epochs} | Train Loss: {train_loss:.4f} Acc: {train_acc:.4f} | Dev Loss: {val_loss:.4f} F1: {val_f1:.4f}")

        if val_f1 > best_dev_f1:
            best_dev_f1 = val_f1
            torch.save(model.state_dict(), best_model_path)

    # Final evaluation on Test set
    if os.path.exists(best_model_path):
        model.load_state_dict(torch.load(best_model_path))

    all_probs, all_labels = [], []
    model.eval()
    with torch.no_grad():
        for sequences, labels in test_loader:
            probs = model(sequences.to(device))
            all_probs.extend(probs.cpu().numpy())
            all_labels.extend(labels.numpy())

    all_probs = np.array(all_probs)
    all_labels = np.array(all_labels)

    np.save(os.path.join(args.output_dir, "lkan_internal_test_probs.npy"), all_probs)
    np.save(os.path.join(args.output_dir, "internal_test_labels.npy"), all_labels)

    preds = (all_probs >= 0.5).astype(int)
    test_acc = accuracy_score(all_labels, preds)
    test_f1 = f1_score(all_labels, preds, zero_division=0)
    test_auc = roc_auc_score(all_labels, all_probs)
    test_prec = precision_score(all_labels, preds, zero_division=0)
    test_rec = recall_score(all_labels, preds, zero_division=0)
    test_mcc = matthews_corrcoef(all_labels, preds)

    print("=" * 50)
    print("FINAL TEST RESULT:")
    print(f"Accuracy: {test_acc:.4f}")
    print(f"F1 Score: {test_f1:.4f} | ROC-AUC: {test_auc:.4f}")
    print(f"Precision: {test_prec:.4f} | Recall: {test_rec:.4f} | MCC: {test_mcc:.4f}")
    print("=" * 50)