# python script to split output of PLINK --recode bimbam by chromosome

# file names from working dir
cord_infile = "cord.recode.geno.txt"
glob_ref_infile = "glob_pure.recode.geno.txt"
glob_mr_infile = "glob_mr.recode.geno.txt"
pos_infile = "glob_mr.recode.pos.txt"
chr_list = ["Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11"]

# loop through each file and write to separate outputs for each chromosome
for chr in chr_list:
    # open input file connections
    cord_in_conn = open(cord_infile, "r")
    glob_ref_in_conn = open(glob_ref_infile, "r")
    glob_mr_in_conn = open(glob_mr_infile, "r")
    pos_in_conn = open(pos_infile, "r")

    # open output file connections
    cord_out_conn = open("geno_by_chr/cord_{CHR}_geno.txt".format(CHR = chr), "w")
    glob_ref_out_conn = open("geno_by_chr/glob_ref_{CHR}_geno.txt".format(CHR = chr), "w")
    glob_mr_out_conn = open("geno_by_chr/glob_mr_{CHR}_geno.txt".format(CHR = chr), "w")
    pos_out_conn = open("geno_by_chr/pos_{CHR}.txt".format(CHR = chr), "w")

    # read through each file, printing headers automatically and otherwise conditionally printing the current loop chromosome
    for line in cord_in_conn:
        if not line.startswith("Chr"):
            cord_out_conn.write(line)
        else:
            if line.startswith(chr):
                cord_out_conn.write(line)
    
    for line in glob_ref_in_conn:
        if not line.startswith("Chr"):
            glob_ref_out_conn.write(line)
        else:
            if line.startswith(chr):
                glob_ref_out_conn.write(line)

    for line in glob_mr_in_conn:
        if not line.startswith("Chr"):
            glob_mr_out_conn.write(line)
        else:
            if line.startswith(chr):
                glob_mr_out_conn.write(line)

    for line in pos_in_conn:
        if line.startswith(chr):
            pos_out_conn.write(line)

    # close input connections
    cord_in_conn.close()
    glob_ref_in_conn.close()
    glob_mr_in_conn.close()
    pos_in_conn.close()

    # close output connections
    cord_out_conn.close()
    glob_ref_out_conn.close()
    glob_mr_out_conn.close()
    pos_out_conn.close()
