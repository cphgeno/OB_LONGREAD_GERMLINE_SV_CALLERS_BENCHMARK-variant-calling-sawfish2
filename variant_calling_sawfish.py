import os
import argparse
import subprocess
import shutil

SEX_INFO_DIR = "./input_data/sawfish_info"

def get_expected_cn_path(sample_name: str) -> str:
    '''
    Read sex information for sample and return expected_cn BED path.
    Looks for .input_data/sawfish_info/<sample_name>.txt

    Arguments:
    ----------
    sample_name: str
        Name of the sample to do the variant calling on.
    '''
    sex_file = os.path.join(SEX_INFO_DIR, f"{sample_name}.txt")

    if not os.path.exists(sex_file):
        raise FileNotFoundError(f"Sex info file not found: {sex_file}")

    with open(sex_file) as f:
        sex = f.read().strip().lower()

    if sex == "female":
        return "/sawfish-v2.2.1-x86_64-unknown-linux-gnu/data/expected_cn/expected_cn.hg38.XX.bed"
    elif sex == "male":
        return "/sawfish-v2.2.1-x86_64-unknown-linux-gnu/data/expected_cn/expected_cn.hg38.XY.bed"
    else:
        raise ValueError(f"Unexpected sex value '{sex}' in {sex_file} (expected 'male' or 'female').")


def variant_calling_sawfish(input_bam: str, reference: str, output_dir: str, name: str) -> None:
    '''
    Performs variant calling with sawfish

    Arguments:
    ----------

    input_bam: str
        input BAM to run sawfish
    output_dir: str
        output storage for vcf
    name: str
        name of analysed sample
    '''
    path_bam = os.path.realpath(input_bam)
    path_discover = os.path.join(output_dir,f"{name}_discover")    
    path_jointcall = os.path.join(output_dir, f"{name}_jointcall")

    reference_real_path = os.path.realpath(reference)

    threads = 16

    # determine expected_cn based on info file
    expected_cn = get_expected_cn_path(name)

    cnv_excluded_regions = "/sawfish-v2.2.1-x86_64-unknown-linux-gnu/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz"

    # sawfish discovery
    discover_cmd = [
        "sawfish",
        "discover",
        "--ref", reference_real_path,
        "--threads", str(threads),
        "--expected-cn", expected_cn,
        "--cnv-excluded-regions", cnv_excluded_regions,
        "--bam", path_bam,
        "--output-dir", path_discover
        ]

    # sawfish joint-call
    joint_call = [
            "sawfish",
            "joint-call",
            "--threads", str(threads),
            "--sample", path_discover,
            "--output-dir", path_jointcall,
            ]

    print("Running:", " ".join(discover_cmd))
    print("Running:", " ".join(joint_call))
    subprocess.run(discover_cmd, check=True)
    subprocess.run(joint_call, check=True)
    
    out_file = f"{output_dir}/{name}.vcf"
    
    # move output to folder
    sawfish_vcf=os.path.join(path_jointcall, "genotyped.sv.vcf.gz") 
    shutil.move(str(sawfish_vcf), str(out_file))
    

def main():
    
    parser = argparse.ArgumentParser(
        description="Run variant calling with Sawfish."
    )

    parser.add_argument(
        "--name",
        type=str,
        required=True,
        help="Sample or dataset name (used as an identifier for outputs)."
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory where output VCFs and intermediate files will be written."
    )

    # Linear alignment inputs
    parser.add_argument(
        "--data.linear.bam",
        dest="data_linear_bam",
        required=True,
        help="Linear-aligned BAM file for the dataset."
    )

    parser.add_argument(
        "--data.linear.bam.bai",
        dest="data_linear_bam_bai",
        required=True,
        help="Index (.bai) for the linear-aligned BAM file."
    )

    # Graph alignment inputs
    parser.add_argument(
        "--data.graph.bam",
        dest="data_graph_bam",
        required=True,
        help="Graph-aligned BAM file for the dataset."
    )

    parser.add_argument(
        "--data.graph.bam.bai",
        dest="data_graph_bam_bai",
        required=True,
        help="Index (.bai) for the graph-aligned BAM file."
    )

    parser.add_argument(
        "--data.gam",
        dest="data_gam",
        required=True,
        help="Graph Alignment Map (GAM) file for the dataset."
    )

    # Assembly inputs
    parser.add_argument(
        "--data.asm.hap1",
        dest="data_asm_hap1",
        required=True,
        help="Haplotype 1 assembly graph (GFA format)."
    )

    parser.add_argument(
        "--data.asm.hap2",
        dest="data_asm_hap2",
        required=True,
        help="Haplotype 2 assembly graph (GFA format)."
    )

    # Truth set
    parser.add_argument(
        "--data.truthset",
        dest="data_truthset",
        required=True,
        help="Truth set VCF (bgzipped) for benchmarking."
    )

    parser.add_argument(
        "--data.truthset.tbi",
        dest="data_truthset_tbi",
        required=True,
        help="Tabix index (.tbi) for the truth set VCF."
    )

    # Reference inputs
    parser.add_argument(
        "--data.reference",
        dest="data_reference",
        required=True,
        help="Reference genome FASTA file."
    )

    parser.add_argument(
        "--data.reference.fai",
        dest="data_reference_fai",
        required=True,
        help="FAI index for the reference genome FASTA."
    )

    # Workflow to follow
    parser.add_argument(
        "--workflow",
        dest="workflow",
        choices=["linear", "graph"],
        required=True,
        help=(
            "Select which BAM alignment to use for variant calling. "
            "'linear' uses --data.linear.bam, "
            "'graph' uses --data.graph.bam."
        )
    )

    args = parser.parse_args()

    if args.workflow == "linear":
        active_bam = args.data_linear_bam
        active_bam_bai = args.data_linear_bam_bai
    elif args.workflow == "graph":
        active_bam = args.data_graph_bam
        active_bam_bai = args.data_graph_bam_bai
    else:
        raise ValueError(f"Unsupported BAM mode: {args.workflow}")

    variant_calling_sawfish(
        active_bam,
        args.data_reference,
        args.output_dir,
        args.name
    )


if __name__ == "__main__":
    main()


