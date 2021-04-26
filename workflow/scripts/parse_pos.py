from parsing.pos_parsing import parse_pos

paths = {
    'bhsadata': snakemake.input.bhsadata,
    'functions': snakemake.input.functions,
    'slot2pos': snakemake.output.slot2pos,
    'uniquepos': snakemake.output.uniquepos
}

# snakemake function execution
parse_pos(paths)
