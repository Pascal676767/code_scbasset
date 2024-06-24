from analyse_function import train_validation_comparison
import argparse

def main():
    parser = argparse.ArgumentParser(description='Compare the training and validation metrics (Loss and AUC).')
    parser.add_argument('file_path', type=str, help='Path to the pickle file containing the metrics data.')
    parser.add_argument('--title_loss', type=str, default='Comparison between train and validation Loss', help='Title for the Loss plot.')
    parser.add_argument('--title_auc', type=str, default='Comparison between train and validation AUC', help='Title for the AUC plot.')
    parser.add_argument('--output_file', type=str, default='./AUC_Loss.png', help='File path to save the output plot.')

    args = parser.parse_args()

    train_validation_comparison(file_path=args.file_path, title_loss=args.title_loss, title_auc=args.title_auc, output_file=args.output_file)

if __name__ == '__main__':
    main()

