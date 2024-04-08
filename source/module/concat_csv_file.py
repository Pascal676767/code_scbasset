def concat_csv_file(file1, file2):
    data_file1 = pd.read_csv(file1, sep = ',', header = 0)
    data_file2 = pd.read_csv(file2, sep = ',', header = 0)
    
    if (data_file1.columns == data_file2.columns).all():
        data_concat = pd.concat([data_file1, data_file2], ignore_index=True)
        return data_concat
        
