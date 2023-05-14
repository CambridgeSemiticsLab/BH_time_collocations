from data.pipeline.legacy_scripts.tools import nav_tree as nt


def get_times(time_parse):
    """Recursively select time words from a parsed timephrase."""
    for time in time_parse['times']:
        if type(time) == int:
            yield time
        elif type(time) == dict:
            yield from get_times(time)
        elif type(time) == list:
            head = nt.get_head(time)
            yield head
        else:
            raise Exception(f'problematic time {time} in {time_parse}')
