# Purpose: Academic-grade CNN-BiLSTM baseline model for DNA sequence classification (e.g., H-DNA prediction).
# Strategy: One-hot encoding + Multi-scale Conv1D (Inception-style) + BiLSTM + Global Max Pooling + BCEWithLogitsLoss.
# Execution Example:
#   python cnn_bilstm_baseline.py --data_dir ./data --output_dir ./output --batch_size 512 --lr 5e-4 --epochs 10

import os
import argparse
from datetime import datetime
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
from sklearn.metrics import average_precision_score, roc_auc_score, matthews_corrcoef, confusion_matrix
from tqdm import tqdm
 

class DNAOneHotDataset(Dataset):
    """
    Dataset wrapper for DNA sequences converting standard ACGT base characters into One-Hot encodings.
    """
    def __init__(self, csv_path, max_len=150):
        df = pd.read_csv(csv_path)
        self.max_len = max_len
        self.sequences = df['sequence'].str.upper().values
        self.labels = df['label'].values
        self.mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        seq = self.sequences[idx][:self.max_len]
        one_hot = np.zeros((4, self.max_len), dtype=np.float32)
        for i, base in enumerate(seq):
            if base in self.mapping:
                one_hot[self.mapping[base], i] = 1.0
        return torch.FloatTensor(one_hot), torch.FloatTensor([self.labels[idx]])


class Academic_CNN_BiLSTM(nn.Module):
    """
    Academic baseline architecture integrating multi-scale 1D convolutions (Inception motif capture),
    a bidirectional LSTM for sequence dependency learning, and global max pooling.
    """
    def __init__(self, hidden_dim=128):
        super().__init__()
        # Multi-scale 1D Convolutional layers to capture motifs of varying lengths
        self.conv7 = nn.Conv1d(4, 64, kernel_size=7, padding=3)
        self.conv9 = nn.Conv1d(4, 64, kernel_size=9, padding=4)
        self.conv11 = nn.Conv1d(4, 64, kernel_size=11, padding=5)

        self.bn1 = nn.BatchNorm1d(192)
        self.relu = nn.ReLU()
        self.pool = nn.MaxPool1d(4)

        # Single-layer Bidirectional LSTM
        self.lstm = nn.LSTM(192, hidden_dim, bidirectional=True, batch_first=True, num_layers=1)

        # Classification head outputting unscaled logits for numeric stability with BCEWithLogitsLoss
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim * 2, 64),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(64, 1)
        )

        self._init_weights()

    def _init_weights(self):
        """Orthogonal initialization for LSTM and Kaiming initialization for Conv layers to ensure robust convergence."""
        for name, param in self.named_parameters():
            if 'conv' in name and 'weight' in name:
                nn.init.kaiming_normal_(param)
            elif 'lstm' in name and 'weight' in name:
                nn.init.orthogonal_(param)
            elif 'bias' in name:
                nn.init.constant_(param, 0)

    def forward(self, x):
        # 1. Multi-scale feature extraction
        x7 = self.conv7(x)
        x9 = self.conv9(x)
        x11 = self.conv11(x)

        # 2. Concatenation and max pooling
        x = torch.cat([x7, x9, x11], dim=1)
        x = self.pool(self.relu(self.bn1(x)))

        # 3. Permute dimensions for LSTM (Batch, Seq_Len, Channel)
        x = x.transpose(1, 2)
        x, _ = self.lstm(x)

        # 4. Global Max Pooling across sequence dimension to capture prominent motif signals
        x = torch.max(x, dim=1)[0]
        return self.classifier(x)


def parse_args():
    """Parse command line arguments for flexible GitHub execution."""
    parser = argparse.ArgumentParser(description="Academic CNN-BiLSTM Baseline for DNA Classification")
    parser.add_argument("--data_dir", type=str, default="./data", help="Directory path containing train.csv, dev.csv, and test.csv")
    parser.add_argument("--output_dir", type=str, default="./output", help="Directory path to save model checkpoints and logs")
    parser.add_argument("--max_len", type=int, default=150, help="Maximum sequence length padding/truncation size")
    parser.add_argument("--batch_size", type=int, default=512, help="Batch size for training and evaluation")
    parser.add_argument("--lr", type=float, default=5e-4, help="Initial learning rate for AdamW optimizer")
    parser.add_argument("--epochs", type=int, default=10, help="Total number of training epochs")
    return parser.parse_args()


def train(args):
    """Main training, evaluation, and logging pipeline."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    timestamp = datetime.now().strftime("%m%d_%H%M")
    save_dir = os.path.join(args.output_dir, f"CNN_BiLSTM_Exp_{timestamp}")
    os.makedirs(save_dir, exist_ok=True)

    # Initialize data loaders
    train_ds = DNAOneHotDataset(os.path.join(args.data_dir, "train.csv"), max_len=args.max_len)
    dev_ds = DNAOneHotDataset(os.path.join(args.data_dir, "dev.csv"), max_len=args.max_len)
    test_ds = DNAOneHotDataset(os.path.join(args.data_dir, "test.csv"), max_len=args.max_len)

    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, num_workers=4, pin_memory=True)
    dev_loader = DataLoader(dev_ds, batch_size=args.batch_size, num_workers=4)
    test_loader = DataLoader(test_ds, batch_size=args.batch_size, num_workers=4)

    # Model initialization
    model = Academic_CNN_BiLSTM().to(device)
    optimizer = optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='max', factor=0.5, patience=1)
    criterion = nn.BCEWithLogitsLoss()

    best_auprc = 0.0
    history = []

    # Training loop
    for epoch in range(args.epochs):
        model.train()
        running_loss = 0.0
        pbar = tqdm(train_loader, desc=f"Epoch {epoch + 1}/{args.epochs}", unit="batch")

        for seqs, labels in pbar:
            seqs, labels = seqs.to(device), labels.to(device)
            optimizer.zero_grad()

            logits = model(seqs)
            loss = criterion(logits, labels)
            loss.backward()

            # Gradient clipping to prevent exploding gradients
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

            running_loss += loss.item()
            with torch.no_grad():
                preds = (torch.sigmoid(logits) >= 0.5).float()
                batch_acc = (preds == labels).sum().item() / labels.size(0)

            pbar.set_postfix({
                "Inst_Loss": f"{loss.item():.4f}",
                "Inst_Acc": f"{batch_acc:.2%}",
                "LR": f"{optimizer.param_groups[0]['lr']:.1e}"
            })

        # Validation stage
        model.eval()
        dev_probs, dev_labels = [], []
        with torch.no_grad():
            for seqs, labels in dev_loader:
                logits = model(seqs.to(device))
                probs = torch.sigmoid(logits)
                dev_probs.extend(probs.detach().cpu().numpy())
                dev_labels.extend(labels.numpy())

        cur_auprc = average_precision_score(dev_labels, dev_probs)
        scheduler.step(cur_auprc)

        # Save optimal model weights based on Dev AUPRC performance
        if cur_auprc > best_auprc:
            best_auprc = cur_auprc
            torch.save(model.state_dict(), os.path.join(save_dir, "best_model.pth"))

            # Evaluate on Test set
            test_probs, test_labels = [], []
            with torch.no_grad():
                for seqs, labels in test_loader:
                    logits = model(seqs.to(device))
                    probs = torch.sigmoid(logits)
                    test_probs.extend(probs.detach().cpu().numpy())
                    test_labels.extend(labels.numpy())

            final_auprc = average_precision_score(test_labels, test_probs)
            final_auroc = roc_auc_score(test_labels, test_probs)
            test_preds = (np.array(test_probs) >= 0.5).astype(int)
            mcc = matthews_corrcoef(test_labels, test_preds)
            tn, fp, fn, tp = confusion_matrix(test_labels, test_preds).ravel()
            spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0

            result_snapshot = {
                "epoch": epoch + 1, "AUPRC": final_auprc, "AUROC": final_auroc,
                "MCC": mcc, "Spec": spec, "TP": int(tp), "FP": int(fp)
            }
            history.append(result_snapshot)

            # Log experiment metrics
            with open(os.path.join(save_dir, "experiment_log.txt"), "a", encoding='utf-8') as f:
                f.write(f"[{datetime.now()}] Epoch {epoch + 1}: Best AUPRC={final_auprc:.4f}, AUROC={final_auroc:.4f}, MCC={mcc:.4f}\n")

    # Export metrics history
    if history:
        pd.DataFrame(history).to_csv(os.path.join(save_dir, "training_history.csv"), index=False)


if __name__ == "__main__":
    args = parse_args()
    train(args)