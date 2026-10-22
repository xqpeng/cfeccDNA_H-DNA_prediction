# Introduction
Models for H-DNA Sequence Identification
This fold contains the official PyTorch & Scikit-Learn implementations for benchmarking deep learning architectures and traditional machine learning baselines on high-throughput H-DNA sequence prediction datasets. 
## Repository Structure
The project is organized into modular directories corresponding to each model family, ensuring native compatibility with Hugging Face Transformers, PyTorch, and Scikit-Learn ecosystems:

 .
 ├── CNN_BiLSTM/
 │   ├── cnn_bilstm.py             # Multi-scale Conv1D + BiLSTM baseline script
 │   ├── best_model.pth            # Trained model weights checkpoint
 ├── DNABERT2/
 │   ├── train_DNABERT2.py         # Fine-tuning entry point using HF Trainer & PEFT
 │   ├── config.json               # Model configuration file
 │   ├── tokenizer_config.json     # Tokenizer configuration specs
 ├── LKAN/
 │   ├── KANLinear.py              # PyTorch B-spline KAN layer implementation
 │   ├── run_LKAN.py               # LKAN architecture execution & evaluation pipeline
 │   ├── best_kan_baseline.pth     # Optimal LKAN model checkpoint
 ├── machine_learning.py           # Scalable ML benchmark (RF, XGBoost, LR) with 4-mer features
 ├── requirements.txt              # Environment dependencies
 └── README.md                     # Repository documentation
 
 

## Environment & Installation
Ensure you have Python 3.9+ and CUDA 11.8/12.1 installed. Install all necessary dependencies via:
git clone https://github.com/xqpeng/cfeccDNA_H-DNA_prediction.git
cd cfeccDNA_H-DNA_prediction/Triplex_prediction_models
pip install -r requirements.txt

## Data Preparation
Place your train, validation, and test datasets in the ./data directory. All CSV files must contain at least two columns:
• sequence: Raw DNA nucleotide sequence (e.g., AGCT...)
• label: Binary target indicator (0 or 1)
Expected directory structure:


./data/
├── train.csv
├── dev.csv
└── test.csv

## Usage & Quick Start
1. **Fine-tuning DNABERT-2**
python DNABERT2/train_DNABERT2.py \
    --model_name_or_path zhihan1996/DNABERT-2-117M \
    --data_path ./data \
    --kmer -1 \
    --run_name dnabert2_hdna_finetuned \
    --model_max_length 150 \
    --per_device_train_batch_size 64 \
    --per_device_eval_batch_size 64 \
    --learning_rate 3e-5 \
    --num_train_epochs 5 \
    --bf16 True \
    --dataloader_num_workers 8 \
    --save_steps 2500 \
    --eval_steps 2500 \
    --output_dir output/dnabert2_hdna_finetuned \
    --evaluation_strategy steps \
    --logging_steps 100 \
    --overwrite_output_dir True \
    --remove_unused_columns True \
    --dataloader_pin_memory True \
    --save_total_limit 2


2. **Multi-scale CNN-BiLSTM Baseline**
python CNN_BiLSTM/cnn_bilstm.py \
    --data_dir ./data \
    --output_dir ./output \
    --batch_size 512 \
    --lr 5e-4 \
    --epochs 10 \
    --max_len 150


3. **Kolmogorov-Arnold Network (LKAN) Baseline**
python LKAN/run_LKAN.py \
    --data_dir ./data \
    --output_dir ./output \
    --target_seq_len 150 \
    --batch_size 512 \
    --lr 5e-4 \
    --epochs 10

4. **Classical Machine Learning Benchmark**
python machine_learning.py \
    --train_path ./data/train.csv \
    --test_path ./data/test.csv \
    --output_dir ./output \
    --kmer 4

5. **Evaluation Metrics**
All models evaluate classification performance comprehensively on test sets using the following metrics:
**AUPRC**: Area Under the Precision-Recall Curve (Primary Target Metric)
**AUROC**: Area Under the Receiver Operating Characteristic
**MCC**: Matthews Correlation Coefficient
**Specificity (Spec)**: True Negative Rate
**F1-Score & Accuracy**: Precision, Recall, and Accuracy
Prediction probabilities are automatically dumped to .npz and .npy formats to facilitate rendering publication-grade ROC and PR curves.


