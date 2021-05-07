# modify default BHSA domain features

def permissive_q(node, api):
    """Map domains to be more permissive with Q
    
    We allow narration to be tagged as Q if it is
    embedded in Q. This avoids some overly structuralist
    assumptions about verb forms which insists that, e.g.,
    a wayyiqtol "has to be" narration.
    
    This feature will return either Q or the last value
    of the feature txt on a clause node.
    """
    txt_type = api.F.txt.v(node)
    if 'Q' in txt_type:
        return 'Q'
    else:
        return txt_type[-1]


