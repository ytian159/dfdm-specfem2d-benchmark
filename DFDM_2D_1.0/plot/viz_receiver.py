import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import sys

def main(directory):
    # Replace 'receiver*.out' with the actual file pattern
    file_pattern = os.path.join(directory, 'receiver*.out')

    # Get a list of all files matching the pattern
    file_list = glob.glob(file_pattern)

    if not file_list:
        print(f"No files matching the pattern 'receiver*.out' found in directory: {directory}")
        return

    # Create subplots
    num_files = len(file_list)
    fig, axes = plt.subplots(num_files, 1, figsize=(10, 6 * num_files), squeeze=False)

    # Loop through each file and generate a plot
    for idx, file_name in enumerate(file_list):
        # Read the .out file using pandas
        # Assuming the file is CSV formatted and columns are separated by whitespace or commas
        data = pd.read_csv(file_name, header=None)

        # Check the first few rows to ensure the file is read correctly
        print(f"Preview of {file_name}:")
        print(data.head())

        # Assign column names if the file does not have a header
        # Assuming the file has at least three columns
        data.columns = ['Time', 'Umid1', 'Umid2']

        # Calculate the average of Umid1 and Umid2
        data['Average'] = (data['Umid1'] + data['Umid2']) / 2

        # Plot the average on the corresponding subplot
        ax = axes[idx, 0]
        ax.plot(data['Time'], data['Average'], label='Average (Umid1 & Umid2)', marker='o')

        # Add labels, title, and legend
        ax.set_xlabel('Time')
        ax.set_ylabel('U')
        ax.set_title(f'Receiver Data {file_name}')
        ax.legend()
        ax.grid(True)

    # Adjust layout and show all plots together
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python viz_receiver.py <directory_path>")
    else:
        main(sys.argv[1])