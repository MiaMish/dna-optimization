from datetime import datetime
import click
from codon_remover import *
from splicing import *


@click.command()
@click.option('--dna-seq', help='DNA sequence', prompt='DNA seq (leave empty for random DNA generation)', default="")
@click.option('--optimization', type=click.Choice(['gggg_removal', 'cccc_removal', 'splicing']), prompt='Optimization type')
@click.option('--start-codon', type=click.Choice(['0', '1', '2']), prompt='Start codon', default='0')
@click.option('--output-file', prompt='Output file path', default='out.txt')
@click.option('--seq-name', prompt='Name for seq', default='')
def cli_main(dna_seq, optimization, start_codon, output_file, seq_name):
    click.echo(click.style('Hello World!', fg='green'))
    start_codon_as_int = int(start_codon)
    if dna_seq == "":
        print("Randomizing DNA seq")
        dna_seq = randomize_dna(start_codon_as_int)
    print(f"Using input DNA seq [{trim_string_for_log(dna_seq)}] (length: {len(dna_seq)}) with start codon {start_codon_as_int}")
    normalized_input = normalize_dna(dna_seq)
    is_valid_input = validate_dna_seq(normalized_input)
    if not is_valid_input:
        return
    if optimization == 'gggg_removal':
        optimized_seq = remove_gggg_codons(normalized_input, start_codon_as_int)
    elif optimization == 'cccc_removal':
        optimized_seq = remove_cccc_codons(normalized_input, start_codon_as_int)
    elif optimization == 'splicing':
        optimized_seq = req_splicing(normalized_input, start_codon_as_int)
    is_optimized_valid = validate_optimized_seq(optimized_seq, normalized_input, start_codon_as_int)
    if not is_optimized_valid:
        return
    analyze_optimized_seq(optimized_seq, normalized_input, start_codon_as_int)
    if seq_name == "":
        seq_name = f"DNA seq {datetime.today().strftime('%Y-%m-%d-%H:%M:%S')}"
        print(f"Seq name was not provided, using auto generated name: {seq_name}")
    print(f"Writing results to {output_file}")
    with open(output_file, 'a') as fout:
        fout.write(f"> Original {seq_name}\n")
        fout.write(f"{dna_seq}\n\n")
        fout.write(f"> Optimized {seq_name}\n")
        fout.write(f"{optimized_seq}\n\n")
    return


cli_main()
