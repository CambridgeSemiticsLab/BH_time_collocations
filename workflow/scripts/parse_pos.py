from parsing.pos_parsing import parse_pos

# snakemake function execution
parse_pos(
    snakemake.input.bhsadata,
    snakemake.output.slot2pos,
    snakemake.output.uniquepos
)
