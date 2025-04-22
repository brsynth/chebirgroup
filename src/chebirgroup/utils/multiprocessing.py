from itertools import islice

def chunk_iterable(iterable, chunk_size: int):
    it = iter(iterable)
    return iter(lambda: list(islice(it, chunk_size)), [])
