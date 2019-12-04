# Modules for working with linguistic data

def is_disjoint(ph, tf):
    """Isolate phrases with gaps."""
    L = tf.api.L
    ph = L.d(ph,'word')
    for w in ph:
        if ph[-1] == w:
            break
        elif (ph[ph.index(w)+1] - w) > 1:
            return True
