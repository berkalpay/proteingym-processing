# # INPUT
# 1. .a2m file (updated in this stage)
# 2. align_config file (updated in this stage)


def uppercase_columns(aln):
    """
    (Reverse the lowercase_column function)
    Change a subset of columns to lowercase character
    and replace "." gaps with "-" gaps
    Change all lower case characters to upper case characters

    Parameters
    ----------
    aln : Alignment object

    Returns
    -------
    Alignment
        Alignment with all uppercase columns
        and only "-" gaps
    """
    columns = np.array([s.islower() for s in aln.matrix[0, :]])
    return aln.apply(columns=columns, func=np.char.upper).replace(
        ".", "-", columns=columns
    )


def update_align_cfg(colcov_70_cfg, colcov_0_cfg, colcov_0_alignment_file, seq_len):
    """
    Update the alignment output config file
    from minimum column coverage = 70
    to look more like that from minimum column coverage = 0
    as this is needed for couplings input at minimum column coverage = 0

    Note: This new alignment output config file has to be in the ${prefix}_colcov_0/align directory

    Parameters
    ----------
    colcov_70_cfg (str): align stage output config file name from minimum column coverage = 70
    colcov_0_cfg (str): align stage output config file name from minimum column coverage = 0
    colcov_0_alignment_file (str): file name of hand-made align stage output alignment
        (converted from colcov_70 to colcov_0)
    seq_len (int): full length of the query sequence from the alignment file

    Returns
    -------
    None
    """
    config = read_config_file(colcov_70_cfg, preserve_order=True)
    config["alignment_file"] = colcov_0_alignment_file
    config["num_sites"] = seq_len
    config["segments"][0][5] = list(range(1, seq_len + 1))

    write_config_file(colcov_0_cfg, config)


# test for couplings only
job_name_prefix = (
    f"{output_folder}/{protein}/{protein}_bit_{bitscore}_theta_{theta}_colcov_70"
)
job_name_prefix_colcov0 = (
    f"{output_folder}/{protein}/{protein}_bit_{bitscore}_theta_{theta}_colcov_0"
)
aln_prefix = (
    job_name_prefix + "/align/" + protein + "_bit_" + bitscore + "_theta_" + theta
)
aln_prefix_colcov0 = (
    job_name_prefix_colcov0
    + "/align/"
    + protein
    + "_bit_"
    + bitscore
    + "_theta_"
    + theta
)
colcov70_alignment = aln_prefix + "_colcov_70" + ".a2m"
colcov0_alignment = aln_prefix + "_colcov_0" + ".a2m"
colcov70_alignment_cfg = aln_prefix + "_colcov_70" + "_align.outcfg"
colcov0_alignment_cfg = aln_prefix_colcov0 + "_colcov_0" + "_align.outcfg"
with open(colcov70_alignment, "r") as infile, open(colcov0_alignment, "w") as outfile:
    aln = Alignment.from_file(infile, format="fasta")
    seq_len = aln.L
    out_aln = uppercase_columns(aln)
    out_aln.write(outfile, format="fasta", width=80)
update_align_cfg(
    colcov70_alignment_cfg, colcov0_alignment_cfg, colcov0_alignment, seq_len
)
os.system(
    (
        f"evcouplings --protein {protein} --stages couplings"
        f"-a {colcov0_alignment} -s {fasta_seq_file} --theta {theta} --colcov 0"
        f"--prefix {job_name_prefix_colcov0} configs/couplings.txt --account marks --yolo"
    )
)

# # OUTPUT
# 1. couplings model with score for each position in the query sequence, regardless of upper/lower case
