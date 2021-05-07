def map_book_collections(label_and_ranges, api):
    """Apply a label to a range of books.
    
    Args:
        label_and_ranges: list of 3-tuples consisting of
            (label, beginning book, ending book) where
            beginning/ending books refers to boundaries
            of the range
        api: an instance of Text-Fabric with BHSA loaded

    Returns:
        dict of book name to new book name strs
    """
    label_dict = {}
    T = api.T
    for label, book_start, book_end in label_and_ranges:
        bs_node = T.nodeFromSection((book_start,))
        be_node = T.nodeFromSection((book_end,))
        in_between = list(range(bs_node, be_node+1))
        whole_section = [bs_node] + in_between + [be_node]
        for book_node in whole_section:
            label_dict[T.sectionFromNode(book_node)[0]] = label 
    return label_dict

def get_book_maps(api):
    """Run book maps using a supplied TF api"""
    
    return dict(
        period=map_book_collections([
            ('SBH', 'Genesis', '2_Kings'),
            ('LBH', 'Esther', '2_Chronicles'),
        ], api),
        tripart=map_book_collections([
            ('Law', 'Genesis', 'Deuteronomy'), 
            ('Prophets', 'Joshua', 'Malachi'), 
            ('Writings', 'Psalms', '2_Chronicles')
        ], api),
        super=map_book_collections([
            ('Samuel', '1_Samuel', '2_Samuel'),
            ('Kings', '1_Kings', '2_Kings'),
            ('Chronicles', '1_Chronicles', '2_Chronicles'),
            ('Ezra-Neh', 'Ezra', 'Nehemiah'),
            ('Twelve', 'Hosea', 'Malachi'),
            ('Megilloth', 'Ruth', 'Esther')
        ], api)
    )

etcbc2sbl = {
    'Genesis':'Gen',
    'Exodus':'Exod',
    'Leviticus':'Lev',
    'Numbers':'Num',
    'Deuteronomy':'Deut',
    'Joshua':'Josh',
    'Judges':'Judg',
    '1_Samuel':'1Sam',
    '2_Samuel':'2Sam',
    '1_Kings':'1Kgs',
    '2_Kings':'2Kgs',
    'Isaiah':'Isa',
    'Jeremiah':'Jer',
    'Ezekiel':'Ezek',
    'Hosea':'Hos',
    'Joel':'Joel',
    'Amos':'Amos',
    'Obadiah':'Obad',
    'Jonah':'Jonah',
    'Micah':'Mic',
    'Nahum':'Nah',
    'Habakkuk':'Hab',
    'Zephaniah':'Zeph',
    'Haggai':'Hag',
    'Zechariah':'Zech',
    'Malachi':'Mal',
    'Psalms':'Ps',
    'Job':'Job',
    'Proverbs':'Prov',
    'Ruth':'Ruth',
    'Song_of_songs':'Song',
    'Ecclesiastes':'Eccl',
    'Lamentations':'Lam',
    'Esther':'Esth',
    'Daniel':'Dan',
    'Ezra':'Ezra',
    'Nehemiah':'Neh',
    '1_Chronicles':'1Chr',
    '2_Chronicles':'2Chr'
}

etcbc2abbr = {
    'Genesis':'Gen.',
    'Exodus':'Exod.',
    'Leviticus':'Lev.',
    'Numbers':'Num.',
    'Deuteronomy':'Deut.',
    'Joshua':'Josh.',
    'Judges':'Judg.',
    '1_Samuel':'I Sam.',
    '2_Samuel':'II Sam.',
    '1_Kings':'I Kgs.',
    '2_Kings':'II Kgs.',
    'Isaiah':'Isa.',
    'Jeremiah':'Jer.',
    'Ezekiel':'Ezek.',
    'Hosea':'Hos.',
    'Joel':'Joel',
    'Amos':'Amos',
    'Obadiah':'Obad.',
    'Jonah':'Jonah',
    'Micah':'Mic.',
    'Nahum':'Nah.',
    'Habakkuk':'Hab.',
    'Zephaniah':'Zeph.',
    'Haggai':'Hag.',
    'Zechariah':'Zech.',
    'Malachi':'Mal.',
    'Psalms':'Psa.',
    'Job':'Job',
    'Proverbs':'Prov.',
    'Ruth':'Ruth',
    'Song_of_songs':'Song',
    'Ecclesiastes':'Eccl.',
    'Lamentations':'Lam.',
    'Esther':'Esth.',
    'Daniel':'Dan.',
    'Ezra':'Ezra',
    'Nehemiah':'Neh.',
    '1_Chronicles':'I Chr.',
    '2_Chronicles':'II Chr.'
}
