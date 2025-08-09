import os
import subprocess
import argparse
import random
from multiprocessing import Pool
from datetime import datetime

def print_colorful_message(message, color):
    """
    Print a colorful message to the console.

    Args:
        message (str): The message to be printed.
        color (str): The color code to be applied.
                     'red', 'green', 'yellow', 'blue', 'magenta', 'cyan', 'white'
    """
    colors = {
        'red': '\033[91m',
        'green': '\033[92m',
        'yellow': '\033[93m',
        'blue': '\033[94m',
        'magenta': '\033[95m',
        'cyan': '\033[96m',
        'white': '\033[97m',
    }
    end_color = '\033[0m'
    if color not in colors:
        print("Invalid color. Please choose from 'red', 'green', 'yellow', 'blue', 'magenta', 'cyan', 'white'.")
        return
    colored_message = f"{colors[color]}{message}{end_color}"
    print(colored_message)

def process_sample(sample_id, path6_rcr, path7_ku2, path8_bracken, db_ku, itm_path, num_threads, force, se):
    task_complete_file = os.path.join(path7_ku2, f"{sample_id}.task.complete")
    if os.path.exists(task_complete_file) and os.path.exists(os.path.join(path7_ku2, f"{sample_id}.kraken.report.txt")) and not force:
        print(f">>== Sample {sample_id} processing already completed. Skipping...")
        return
    
    print("  ")
    print(f">>>== Processing sample: {sample_id} .")
    if se:
        faq_mr1 = os.path.join(path6_rcr, f"{sample_id}_rcr.fastq.gz")
        kraken_report = os.path.join(path7_ku2, f"{sample_id}.kraken.report.txt")
        kraken_out = os.path.join(path7_ku2, f"{sample_id}.kraken.output.txt")

        # Run Kraken2
        subprocess.run(["kraken2", "--db", db_ku, "--threads", str(num_threads), "--report-minimizer-data", "--report",
                        kraken_report, "--use-names", "--output", kraken_out, faq_mr1])
    else:
        faq_mr1 = os.path.join(path6_rcr, f"{sample_id}_rcr_1.fastq.gz")
        faq_mr2 = os.path.join(path6_rcr, f"{sample_id}_rcr_2.fastq.gz")
        kraken_report = os.path.join(path7_ku2, f"{sample_id}.kraken.report.txt")
        kraken_out = os.path.join(path7_ku2, f"{sample_id}.kraken.output.txt")

        # Run Kraken2
        subprocess.run(["kraken2", "--db", db_ku, "--threads", str(num_threads), "--paired", "--report-minimizer-data", "--report",
                        kraken_report, "--use-names", "--output", kraken_out, "--paired", faq_mr1, faq_mr2])

    # Convert Kraken report to MPA format
    kraken_report_std = os.path.join(path7_ku2, f"{sample_id}.kraken.report.std.txt")
    with open(kraken_report_std, "w") as f:
        subprocess.run(["cut", "-f1-3,6-8", kraken_report], stdout=f)
    mpa_output_std = os.path.join(path7_ku2, f"{sample_id}.kraken.mpa.std.txt")
    with open(mpa_output_std, "w") as f:
        subprocess.run(["python", os.path.join(itm_path, "itm_helper", "kreport2mpa.py"), "-r", kraken_report_std, "-o", mpa_output_std], stdout=f)

    # Run Bracken for abundance estimation
    bracken_output_g = os.path.join(path8_bracken, f"{sample_id}.g.bracken")
    bracken_output_s = os.path.join(path8_bracken, f"{sample_id}.s.bracken")
    bracken_output_f = os.path.join(path8_bracken, f"{sample_id}.f.bracken")
    bracken_output_o = os.path.join(path8_bracken, f"{sample_id}.o.bracken")
    subprocess.run(["bracken", "-d", db_ku, "-i", kraken_report, "-o", bracken_output_g, "-r", "100", "-l", "G", "-t", "2"])
    subprocess.run(["bracken", "-d", db_ku, "-i", kraken_report, "-o", bracken_output_s, "-r", "100", "-l", "S", "-t", "2"])
    subprocess.run(["bracken", "-d", db_ku, "-i", kraken_report, "-o", bracken_output_f, "-r", "100", "-l", "F", "-t", "2"])
    subprocess.run(["bracken", "-d", db_ku, "-i", kraken_report, "-o", bracken_output_o, "-r", "100", "-l", "O", "-t", "2"])

    # Calculate alpha diversity
    diversity_output_g = os.path.join(path8_bracken, f"{sample_id}.diversity.g.txt")
    diversity_output_s = os.path.join(path8_bracken, f"{sample_id}.diversity.s.txt")
    itm_path2 = os.path.join(itm_path, "itm_helper")

    with open(diversity_output_g, "a") as f:
        subprocess.run(["python", os.path.join(itm_path2, "alpha_diversity.py"), "-f", bracken_output_g, "-a", "Sh"], stdout=f)
        subprocess.run(["python", os.path.join(itm_path2, "alpha_diversity.py"), "-f", bracken_output_g, "-a", "BP"], stdout=f)
        subprocess.run(["python", os.path.join(itm_path2, "alpha_diversity.py"), "-f", bracken_output_g, "-a", "Si"], stdout=f)
        subprocess.run(["python", os.path.join(itm_path2, "alpha_diversity.py"), "-f", bracken_output_g, "-a", "ISi"], stdout=f)
        subprocess.run(["python", os.path.join(itm_path2, "alpha_diversity.py"), "-f", bracken_output_g, "-a", "F"], stdout=f)
    with open(diversity_output_s, "a") as f:
        subprocess.run(["python", os.path.join(itm_path2, "alpha_diversity.py"), "-f", bracken_output_s, "-a", "Sh"], stdout=f)
        subprocess.run(["python", os.path.join(itm_path2, "alpha_diversity.py"), "-f", bracken_output_s, "-a", "BP"], stdout=f)
        subprocess.run(["python", os.path.join(itm_path2, "alpha_diversity.py"), "-f", bracken_output_s, "-a", "Si"], stdout=f)
        subprocess.run(["python", os.path.join(itm_path2, "alpha_diversity.py"), "-f", bracken_output_s, "-a", "ISi"], stdout=f)
        subprocess.run(["python", os.path.join(itm_path2, "alpha_diversity.py"), "-f", bracken_output_s, "-a", "F"], stdout=f)

    # Mark sample as completed
    with open(task_complete_file, "w") as f:
        f.write(">>>== Processing completed.")

def step6_kraken2Bracken(path6_rcr, path7_ku2, path8_bracken, db_ku, itm_path, num_threads=8, force=False, se=False):
    print("   ")
    print_colorful_message("#########################################################", "blue")
    print_colorful_message(" ITMfinder: Identifing Intratumoral Microbiome pipeline ", "cyan")
    print_colorful_message(" If you encounter any issues, please report them at ", "cyan")
    print_colorful_message(" https://github.com/LiaoWJLab/ITMfinder/issues ", "cyan")
    print_colorful_message("#########################################################", "blue")
    print(" Author: Dongqiang Zeng, Qianqian Mao ")
    print(" Email: interlaken@smu.edu.cn ")
    print_colorful_message("#########################################################", "blue")
    print("   ")

    print(" >>> Perform taxonomic classification using Kraken2...")
    # Create output directories if they do not exist
    os.makedirs(path7_ku2, exist_ok=True)
    os.makedirs(path8_bracken, exist_ok=True)

    # Get list of microbiome reads files
    if se:
         mr_files = [file for file in os.listdir(path6_rcr) if file.endswith("_rcr.fastq.gz")]
    else:
         mr_files = [file for file in os.listdir(path6_rcr) if file.endswith("_rcr_1.fastq.gz")]

    total_files = len(mr_files)
    random.shuffle(mr_files)  # Shuffle the list of files

    for sample_id in mr_files:
       if se:
           sample_id = sample_id[:-len("_rcr.fastq.gz")]  # Remove "_rcr.fastq.gz" suffix
       else:
          sample_id = sample_id[:-len("_rcr_1.fastq.gz")]  # Remove "_rcr_1.fastq.gz" suffix
    
       process_sample(sample_id, path6_rcr, path7_ku2, path8_bracken, db_ku, itm_path, num_threads, force, se)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 6: Taxonomic classification using Kraken2 for report")
    parser.add_argument("--path6_rcr", type=str, help="Path to microbiome reads after decontainmination")
    parser.add_argument("--path7_ku2", type=str, help="Path to Kraken2 outputs for report")
    parser.add_argument("--path8_bracken", type=str, help="Path for Bracken outputs")
    parser.add_argument("--db_ku", type=str, help="Path to Kraken2 database")
    parser.add_argument("--itm_path", type=str, help="Path to ITMfinder")
    parser.add_argument("--num_threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--force", action="store_true", help="Force execution of downstream Kraken2 process")
    parser.add_argument("--se", action="store_true", help="Indicate single-end reads")
    args = parser.parse_args()

    # Call the function with provided arguments
    step6_kraken2Bracken(args.path6_rcr, args.path7_ku2, args.path8_bracken, args.db_ku, args.itm_path, args.num_threads, args.force, args.se)
