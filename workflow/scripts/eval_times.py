from parsing.time_evaler import time_evaler

paths = { 
    'bhsadata': snakemake.input.bhsadata,
    'corrections': snakemake.input.corrections,
    'parsed': snakemake.input.parsed,
    'translations': snakemake.input.translations,
    'ntoget': snakemake.params.ntoget,
    'todo': snakemake.output.todo,
    'evaled': snakemake.output.evaled,
}

time_evaler(
    paths,
)
