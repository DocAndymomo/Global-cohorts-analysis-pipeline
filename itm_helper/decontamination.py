conta_file = os.path.join(itm_path, "itm_helper","conta_ls.txt")
    
    with open(conta_file, "r") as file:
         lines = file.readlines()

    second_column = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) > 1:
            second_column.append(parts[1])

    tls = " ".join(second_column)
    conta_ls = tls.split()
    
    if is_single_end:
        process = subprocess.Popen(["python", os.path.join(itm_path, "itm_helper", "extract_kraken_reads.py"),
                                    "-k", kraken_output, "-U", output_r1, "-o", output_r3,
                                    "--taxid", *conta_ls, "--exclude", "--include-children", "-r", report_file])

    else:
        process = subprocess.Popen(["python", os.path.join(itm_path, "itm_helper", "extract_kraken_reads.py"),
                                    "-k", kraken_output, "-s1", output_r1, "-s2", output_r2, "-o", output_r3, "-o2", output_r4,
                                    "--taxid", *conta_ls, "--exclude", "--include-children", "-r", report_file])
    
