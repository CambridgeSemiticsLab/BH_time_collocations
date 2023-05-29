"""Configure the build of BHSA KinghamThesis version."""

from kingham_thesis.data_pipeline.corpus_building.utils import EditAction

# metadata unique to annotation features, to override the default etcbc values
annotation_metadata = {
    'encoders': 'Cody Kingham',
    'source': ('Kingham, Constructional Network of Time Adverbials in Biblical Hebrew, '
               'PhD Thesis, Cambridge University'),
    'source-url': 'https://github.com/CambridgeSemiticsLab/BH_time_collocations',
}


# Corpus edits are logged here;
# NB: All node numbers are from ETCBC 2021 version;
# New nodes will be assigned with the following scheme:
#     phrases = 2000000 + N
#     phrase_atoms = 2100000 + N
#     subphrases = 2200000 + N
# Note that these temporary node numbers will be re-assigned when
# the corpus is re-indexed during the build process

THESIS_CORPUS_PARAMS = dict(
    book_limit='2_Kings',
    delete_features={
        'book@am', 'book@ar', 'book@bn', 'book@da',
        'book@de', 'book@el', 'book@es', 'book@fa',
        'book@fr', 'book@he', 'book@hi', 'book@id',
        'book@ja', 'book@ko', 'book@la', 'book@nl',
        'book@pa', 'book@pt', 'book@ru', 'book@sw',
        'book@syc', 'book@tr', 'book@ur', 'book@yo',
        'book@zh', 'dist', 'dist_unit',
        'mother_object_type', 'functional_parent',
        'distributional_parent', 'languageISO',
        'omap@c-KT', 'omap@c-2021', 'omap@2017-2021',
    },
    rename_features={},
    edit_actions=[
        EditAction(
            edge_updates={
                'oslots': {
                    661278: {15859, 15860, 15861, 15862},
                    914944: {15859, 15860, 15861, 15862},
                    2200000: {15861},
                    2200001: {15862},
                },
                'mother': {
                    2200001: {2200000},
                },
            },
            feature_updates={
                'otype': {
                    2200000: 'subphrase',
                    2200001: 'subphrase',
                },
                'rela': {
                    2200000: 'NA',
                    2200001: 'dem',
                },
                'typ': {
                    430754: 'Way0',
                    518984: 'Way0',
                },
            },
            deletions={661279, 914945},
            description="Gen 30:16; BLJLH H>W, merge phrases",
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    662567: {17889, 17890, 17891, 17892},
                    916291: {17889, 17890, 17891, 17892},
                    2200002: {17891},
                    2200003: {17892},
                },
                'mother': {
                    2200003: {2200002},
                },
            },
            feature_updates={
                'otype': {
                    2200002: 'subphrase',
                    2200003: 'subphrase',
                },
                'typ': {
                    431180: 'Way0',
                    519421: 'Way0',
                },
                'rela': {
                    2200002: 'NA',
                    2200003: 'dem',
                },
            },
            deletions={662568, 916292},
            description="Gen 32:23;  BLJLH H>W, merge phrases",
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    # clause
                    435375: {
                        37920, 37921, 37922, 37923, 37924,
                        37925, 37926, 37927, 37928, 37929,
                    },
                    # clause_atom
                    523704: {
                        37920, 37921, 37922, 37923, 37924,
                        37925, 37926, 37927, 37928, 37929,
                    },
                    # sentence
                    1178355: {
                        37920, 37921, 37922, 37923, 37924,
                        37925, 37926, 37927, 37928, 37929,
                    },
                },
                'mother': {
                    929211: {929210},  # phrase_atom
                    523704: {523702},  # clause_atom
                },
            },
            deletions={435374, 523703, 1178354},
            feature_updates={
                'rela': {
                    929211: 'Appo',
                },
                'typ': {
                    435375: 'WxY0',
                    523704: 'WxY0',
                },
                'code': {
                    523704: 411,
                },
                'tab': {
                    523704: 6,
                },
            },
            description='Ex 16:26; delete clause/clause atom & merge of its elements',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    438428: {
                        55605, 55606, 55607, 55608, 55609, 55610
                    },
                    526878: {
                        55605, 55606, 55607, 55608, 55609, 55610
                    },
                    438427: {
                        55604, 55611, 55612, 55613, 55614,
                        55615, 55616, 55617, 55618, 55619
                    },
                    526879: {
                        55611, 55612, 55613, 55614,
                        55615, 55616, 55617, 55618, 55619
                    },
                },
            },
            description='Lev 7:17; redraw clause boundaries to include time phrase'
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    688313: {
                        62875, 62876, 62877, 62878, 62879
                    },
                    # NB: use 2M+N for new nodes;
                    # these large numbers will get reassigned later
                    2000000: {
                        62880, 62881, 62882, 62883, 62884, 62885,
                    },
                },
                'head': {
                    62875: {688313},
                    62880: {2000000},
                },
                'nhead': {
                    62877: {688313},
                    62882: {2000000},
                },
            },
            feature_updates={
                'otype': {2000000: 'phrase'},
                'function': {2000000: 'Time'},
                'typ': {2000000: 'PP'},
            },
            description='Lev 16:29; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    690492: {
                        66514, 66515, 66516, 66517, 66518,
                    },
                    2000001: {
                        66519, 66520, 66521,
                    },
                    2000002: {66522, 66523},
                    2000003: {66524, 66525},
                },
                'head': {
                    66519: {2000001},
                    66522: {2000002},
                    66524: {2000003},
                },
                'nhead': {
                    66521: {2000001},
                    66523: {2000002},
                    66525: {2000003},
                },
            },
            feature_updates={
                'otype': {
                    2000001: 'phrase',
                    2000002: 'phrase',
                    2000003: 'phrase',
                },
                'function': {
                    2000001: 'Time',
                    2000002: 'Time',
                    2000003: 'Time',
                },
                'typ': {
                    2000001: 'PP',
                    2000002: 'PP',
                    2000003: 'PP',
                },
            },
            description='Lev 23:32; split up phrase'
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    690894: {67313, 67314, 67315, 67316, 67317},
                    2000004: {67318, 67319, 67320, 67321, 67322, 67323},
                },
                'head': {67318: {2000004}},
                'nhead': {67320: {2000004}},
            },
            feature_updates={
                'otype': {2000004: 'phrase'},
                'function': {2000004: 'Time'},
                'typ': {2000004: 'PP'},
            },
            description='Lev 25:9; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    692414: {
                        69623, 69624, 69625, 69626,
                        69627, 69628, 69629
                    },
                    2000006: {
                        69630, 69631, 69632, 69633, 69634
                    },
                },
                'head': {69630: {2000006}},
                'nhead': {69632: {2000006}},
            },
            feature_updates={
                'otype': {2000006: 'phrase'},
                'function': {2000006: 'Time'},
                'typ': {2000006: 'PP'},
            },
            description='Num 1:1; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    694945: {
                        75786, 75787, 75788, 75789,
                    },
                    2000007: {
                        75790, 75791, 75792, 75793, 75794,
                    },
                },
                'head': {75790: {2000007}},
                'nhead': {75792: {2000007}},
            },
            feature_updates={
                'otype': {2000007: 'phrase'},
                'function': {2000007: 'Time'},
                'typ': {2000007: 'PP'},
            },
            description='Num 9:3; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    694962: {
                        75826, 75827, 75828
                    },
                    2000008: {
                        75829, 75830, 75831, 75832,
                        75833, 75834, 75835,
                    },
                },
                'head': {75829: {2000008}},
                'nhead': {75832: {2000008}},
            },
            feature_updates={
                'otype': {2000008: 'phrase'},
                'function': {2000008: 'Time'},
                'typ': {2000008: 'PP'},
            },
            description='Num 9:5; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    695033: {
                        75958, 75959, 75960, 75961, 75962,
                    },
                    2000009: {
                        75963, 75964, 75965, 75966
                    },
                },
                'head': {75963: {2000009}},
                'nhead': {75966: {2000009}},
            },
            feature_updates={
                'otype': {2000009: 'phrase'},
                'function': {2000009: 'Time'},
                'typ': {2000009: 'PP'},
            },
            description='Num 9:11; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    695308: {
                        76432, 76433, 76434, 76435, 76436
                    },
                    2000010: {
                        76437, 76438, 76439, 76440, 76441
                    },
                    2000011: {
                        76442, 76443, 76444, 76445, 76446
                    },
                },
                'head': {
                    76437: {2000010},
                    76442: {2000011},
                },
                'nhead': {
                    76439: {2000010},
                    76443: {2000011},
                },
            },
            feature_updates={
                'otype': {
                    2000010: 'phrase',
                    2000011: 'phrase',
                },
                'function': {
                    2000010: 'Time',
                    2000011: 'Time',
                },
                'typ': {
                    2000010: 'PP',
                    2000011: 'PP',
                },
            },
            description='Num 10:11; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    # phrases
                    698727: {
                        82437, 82438, 82439, 82440, 82441
                    },
                    2000012: {82442},
                    2000013: {
                        82443, 82444, 82445, 82446, 82447,
                    },
                    # phrase atoms
                    954435: {
                        82437, 82438, 82439, 82440, 82441,
                    },
                    2100012: {82442},
                    2100013: {
                        82443, 82444, 82445, 82446, 82447,
                    },
                },
                'head': {
                    82437: {698727},
                    82442: {2000012},
                    82443: {2000013},
                },
                'nhead': {
                    82439: {698727},
                    82445: {2000013},
                },
            },
            feature_updates={
                'otype': {
                    2000012: 'phrase',
                    2000013: 'phrase',
                    2100012: 'phrase_atom',
                    2100013: 'phrase_atom',
                },
                'function': {
                    2000012: 'Conj',
                    2000013: 'Time',
                },
                'typ': {
                    2000012: 'CP',
                    2000013: 'PP',
                    2100012: 'CP',
                    2100013: 'PP',
                },
                'rela': {
                    2100012: 'NA',
                    2100013: 'NA',
                },
            },
            description='Num 9:19; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    701823: {
                        87869, 87870, 87871, 87872, 87873,
                    },
                    2000014: {
                        87874, 87875, 87876, 87877, 87878,
                    },
                },
                'head': {87874: {2000014}},
                'nhead': {87875: {2000014}},
            },
            feature_updates={
                'otype': {2000014: 'phrase'},
                'function': {2000014: 'Time'},
                'typ': {2000014: 'PP'},
            },
            description='Num 29:1; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    # phrases
                    702396: {
                        89241, 89242, 89243, 89244, 89245,
                    },
                    2000015: {89246},
                    2000016: {
                        89247, 89248, 89249, 89250, 89251
                    },

                    # phrase atoms
                    958402: {
                        89241, 89242, 89243, 89244, 89245,
                    },
                    2100015: {89246},
                    2100016: {
                        89247, 89248, 89249, 89250, 89251
                    }
                },
                'head': {
                    89241: {702396},
                    89246: {2000015},
                    89247: {2000016},
                },
                'nhead': {
                    89243: {702396},
                    89249: {2000016},
                },
            },
            feature_updates={
                'otype': {
                    2000015: 'phrase',
                    2000016: 'phrase',
                    2100015: 'phrase_atom',
                    2100016: 'phrase_atom',
                },
                'function': {
                    2000015: 'Conj',
                    2000016: 'Time',
                },
                'typ': {
                    2000015: 'CP',
                    2000016: 'PP',
                    2100015: 'CP',
                    2100016: 'PP',
                },
                'rela': {
                    2100015: 'NA',
                    2100016: 'NA',
                },
            },
            description='Num 31:19; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    703048: {
                        90705, 90706, 90707, 90708, 90709,
                    },
                    2000017: {
                        90710, 90711, 90712, 90713, 90714,
                        90715, 90716, 90717, 90718
                    },
                },
                'head': {90710: {2000017}},
                'nhead': {90713: {2000017}},
            },
            feature_updates={
                'otype': {2000017: 'phrase'},
                'function': {2000017: 'Time'},
                'typ': {2000017: 'PP'},
            },
            description='Num 33:3; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    703305: {
                        91100, 91101, 91102, 91103,
                    },
                    2000018: {
                        91111, 91112, 91113, 91114, 91115,
                    },
                    2000019: {
                        91116, 91117, 91118, 91119, 91120,
                    },
                },
                'head': {
                    91111: {2000018},
                    91116: {2000019},
                },
                'nhead': {
                    91113: {2000018},
                    91117: {2000019},
                },
            },
            feature_updates={
                'otype': {
                    2000018: 'phrase',
                    2000019: 'phrase',
                },
                'function': {
                    2000018: 'Time',
                    2000019: 'Time',
                },
                'typ': {
                    2000018: 'PP',
                    2000019: 'PP',
                },
            },
            description='Num 33:38; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    # clause
                    447183: {
                        105113, 105114, 105115, 105116
                    },
                    447184: {
                        105124, 105125,
                    },
                    # clause_atom
                    536021: {
                        105113, 105114, 105115, 105116
                    },
                    # sentence
                    1186375: {
                        105113, 105114, 105115, 105116,
                        105117, 105118, 105119, 105120,
                        105121, 105122, 105123,
                    },
                    1186376: {
                        105124, 105125, 105126, 105127,
                        105128, 105129, 105130, 105131,
                        105132, 105133, 105134, 105135,
                        105136, 105137, 105138,
                    },
                },
                'mother': {
                    536023: {536021},
                    536025: {536021},
                    536027: {536021},
                },
            },
            feature_updates={
                'typ': {
                    447183: 'WQt0',
                    536021: 'WQt0',
                },
            },
            deletions={536022},
            description='Deut 21:16; redraw clause lines and move BJWM to prev clause',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    717027: {
                        115204, 115205, 115206, 115207,
                        115208, 115209, 115210,
                    },
                    2000020: {
                        115211, 115212, 115213,
                    },
                },
                'head': {115211: {2000020}},
                'nhead': {115213: {2000020}},
            },
            feature_updates={
                'otype': {2000020: 'phrase'},
                'function': {2000020: 'Time'},
                'typ': {2000020: 'PP'},
            },
            description='Josh 5:10; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    # clause
                    452315: {
                        133556, 133557, 133558, 133559, 133560
                    },
                    452316: {
                        133565,
                    },
                    # clause_atom
                    541303: {
                        133556, 133557, 133558, 133559, 133560
                    },
                    # sentence
                    1189949: {
                        133556, 133557, 133558, 133559, 133560,
                        133561, 133562, 133563, 133564
                    },
                    1189950: {
                        133565
                    },
                },
                'mother': {
                    541305: {541303},
                    541306: {541303},
                    541307: {541303},
                },
            },
            feature_updates={
                'typ': {
                    452315: 'WQt0',
                    541303: 'WQt0',
                    452316: 'ZYq0',
                },
            },
            deletions={541304},
            description='Judg 9:33; redraw clause bounds to include BBQR with WHJH',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    730312: {
                        139255, 139256, 139257
                    },
                    2000021: {
                        139258, 139259, 139260, 139261, 139262
                    },
                },
                'head': {
                    139258: {2000021},
                },
                'nhead': {
                    139260: {2000021},
                },
            },
            feature_updates={
                'otype': {2000021: 'phrase'},
                'function': {2000021: 'Time'},
                'typ': {2000021: 'PP'},
            },
            description='Judg 19:8; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    739001: {
                        153077, 153078, 153079, 153080,
                    },
                    997288: {
                        153077, 153078, 153079, 153080,
                    },
                    2200004: {153079},
                    2200005: {153080},
                },
                'mother': {
                   2200005: {2200004},
                },
            },
            deletions={739002, 997289},
            feature_updates={
                'otype': {
                    2200004: 'subphrase',
                    2200005: 'subphrase',
                },
                'rela': {
                    2200004: 'NA',
                    2200005: 'dem',
                },
                'typ': {
                    456392: 'Way0',
                    545456: 'Way0',
                },
            },
            description='1 Sam 19:10; merge two phrases, keeping HW> as demonstrative'
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    747140: {
                        166007, 166008, 166009, 166010, 166011
                    },
                    2000022: {
                        166012,
                    },
                    2000023: {
                        166013, 166014,
                    },
                },
                'head': {
                    166007: {747140},
                    166012: {2000022},
                    166013: {2000023},
                },
                'nhead': {
                    166009: {747140},
                    166014: {2000023},
                },
            },
            feature_updates={
                'otype': {
                    2000022: 'phrase',
                    2000023: 'phrase',
                },
                'function': {
                    2000022: 'Conj',
                    2000023: 'Time',
                },
                'typ': {
                    2000022: 'CP',
                    2000023: 'PP',
                },
            },
            description='1 Sam 11:12; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    752165: {
                        173952, 173953, 173954,
                    },
                    2000024: {
                        173955, 173956, 173957,
                    },
                    2000025: {
                        173958, 173959, 173960, 173961,
                    },
                },
                'head': {
                    173955: {2000024},
                    173958: {2000025},
                },
                'nhead': {
                    173957: {2000024},
                    173959: {2000025},
                },
            },
            feature_updates={
                'otype': {
                    2000024: 'phrase',
                    2000025: 'phrase',
                },
                'function': {
                    2000024: 'Time',
                    2000025: 'Time',
                },
                'typ': {
                    2000024: 'PP',
                    2000025: 'PP',
                },
            },
            description='2 Sam 21:9; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    # clause
                    461570: {
                        178174, 178175, 178176, 178177,
                    },
                    461571: {
                        178184, 178185,
                    },
                    # clause_atom
                    550740: {
                        178174, 178175, 178176, 178177,
                    },
                    # sentence
                    1197110: {
                        178174, 178175, 178176, 178177,
                        178178, 178179, 178180, 178181,
                        178182, 178183,
                    },
                    1197111: {
                        178184, 178185, 178186, 178187, 178188
                    },
                },
                'mother': {
                    550742: {550740},
                    550744: {550740},
                },
            },
            feature_updates={
                'typ': {
                    461570: 'WQt0',
                    550740: 'WQt0',
                },
            },
            deletions={550741},
            description='1 Kgs 2:37; redraw clause lines to include BJWM with WHJH',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    755711: {
                        180079, 180080, 180081, 180082, 180083,
                    },
                    2000026: {
                        180084, 180085, 180086,
                    },
                },
                'head': {
                    180084: {2000025},
                },
                'nhead': {
                    180085: {2000025},
                },
            },
            feature_updates={
                'otype': {2000026: 'phrase'},
                'function': {2000026: 'Time'},
                'typ': {2000026: 'PP'},
            },
            description='1 Kgs 6:1; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    756053: {
                        180819, 180820, 180821, 180822, 180823, 180824,
                    },
                    2000027: {
                        180825, 180826, 180827,
                    },
                },
                'head': {
                    180825: {2000027},
                },
                'nhead': {
                    180826: {2000027},
                },
            },
            feature_updates={
                'otype': {2000027: 'phrase'},
                'function': {2000027: 'Time'},
                'typ': {2000027: 'PP'},
            },
            description='1 Kgs 6:38; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    756573: {
                        182107, 182108, 182109, 182110,
                    },
                    2000028: {
                        182111, 182112, 182113,
                    },
                },
                'head': {
                    182111: {2000028},
                },
                'nhead': {
                    182113: {2000028},
                },
            },
            feature_updates={
                'otype': {2000028: 'phrase'},
                'function': {2000028: 'Time'},
                'typ': {2000028: 'PP'},
            },
            description='1 Kgs 8:2; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    759061: {
                        186694, 186695, 186696, 186697, 186698
                    },
                    2000029: {
                        186699, 186700, 186701, 186702, 186703,
                        186704, 186705,
                    },
                },
                'head': {
                    186699: {2000029},
                },
                'nhead': {
                    186702: {2000029},
                },
            },
            feature_updates={
                'otype': {2000029: 'phrase'},
                'function': {2000029: 'Time'},
                'typ': {2000029: 'PP'},
            },
            description='1 Kgs 12:32; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    759087: {
                        186747, 186748, 186749, 186750,
                    },
                    2000030: {
                        186751, 186752, 186753, 186754, 186755,
                    },
                    2000031: {
                        186756, 186757, 186758,
                    },
                },
                'head': {
                    186751: {2000030},
                    186756: {2000031},
                },
                'nhead': {
                    186753: {2000030},
                    186758: {2000031},
                },
            },
            feature_updates={
                'otype': {
                    2000030: 'phrase',
                    2000031: 'phrase',
                },
                'function': {
                    2000030: 'Time',
                    2000031: 'Time',
                },
                'typ': {
                    2000030: 'PP',
                    2000031: 'PP',
                },
            },
            description='1 Kgs 12:33; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    # phrase
                    773807: {
                        211330, 211331, 211332, 211333, 211334,
                    },
                    2000032: {
                        211335, 211336, 211337, 211338, 211339,
                        211340,
                    },
                    # clause
                    468054: {
                        211322, 211323, 211324, 211325, 211326,
                        211327, 211330, 211331, 211332, 211333,
                        211334, 211335, 211336, 211337, 211338,
                        211339, 211340,
                    },
                    468056: {
                        211341, 211342, 211343, 211344, 211345,
                        211346, 211347, 211348, 211349, 211350,
                    },
                    # clause_atom
                    557416: {
                        211322, 211323, 211324, 211325, 211326,
                        211327, 211330, 211331, 211332, 211333,
                        211334, 211335, 211336, 211337, 211338,
                        211339, 211340,
                    },
                    557418: {
                        211341, 211342, 211343, 211344, 211345,
                        211346, 211347, 211348, 211349, 211350,
                    },
                    # sentence
                    1201759: {
                        211322, 211323, 211324, 211325, 211326,
                        211327, 211328, 211329, 211330, 211331,
                        211332, 211333, 211334, 211335, 211336,
                        211337, 211338, 211339, 211340,
                    },
                    1201760: {
                        211341, 211342, 211343, 211344, 211345,
                        211346, 211347, 211348, 211349, 211350,
                    },
                },
                'head': {
                    211335: {2000032},
                },
                'nhead': {
                    211337: {2000032},
                },
            },
            feature_updates={
                'otype': {
                    2000032: 'phrase',
                },
                'function': {
                    2000032: 'Time',
                },
                'typ': {
                    468056: 'ZQtX',
                    557418: 'ZQtX',
                    2000032: 'PP',
                },
            },
            description='2 Kgs 25:1; redraw clause lines and split phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    773888: {
                        211482, 211483, 211484, 211485, 211486,
                    },
                    2000033: {
                        211487, 211488, 211489, 211490, 211491,
                    },
                },
                'head': {
                    211487: {2000033},
                },
                'nhead': {
                    211488: {2000033},
                },
            },
            feature_updates={
                'otype': {2000033: 'phrase'},
                'function': {2000033: 'Time'},
                'typ': {2000033: 'PP'},
            },
            description='1 Kgs 25:8; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    774087: {
                        211989, 211990, 211991, 211992, 211993, 211994,
                        211995, 211996, 211997, 211998
                    },
                    2000034: {
                        211999, 212000, 212001, 212002,
                    },
                    2000035: {
                        212003, 212004, 212005, 212006, 212007, 212008,
                        212009
                    },
                },
                'head': {
                    211999: {2000034},
                    212003: {2000035},
                },
                'nhead': {
                    212002: {2000034},
                    212004: {2000035},
                },
            },
            feature_updates={
                'otype': {
                    2000034: 'phrase',
                    2000035: 'phrase',
                },
                'function': {
                    2000034: 'Time',
                    2000035: 'Time',
                },
                'typ': {
                    2000034: 'PP',
                    2000035: 'PP',
                },
            },
            description='2 Kgs 25:27; split up phrase',
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    # phrases
                    680921: {
                        48203, 48204, 48205
                    },
                    2000036: {
                        48206,
                    },
                    2000037: {
                        48207, 48208, 48209
                    },
                    # phrase_atoms
                    935296: {
                        48203, 48204, 48205
                    },
                    2100036: {
                        48206,
                    },
                    2100037: {
                        48207, 48208, 48209
                    },
                },
                'head': {
                    48203: {680921},
                    48206: {2000036},
                    48207: {2000037},
                },
                'nhead': {
                    48205: {680921},
                    48209: {2000037},
                },
            },
            feature_updates={
                'otype': {
                    2000036: 'phrase',
                    2000037: 'phrase',
                    2100036: 'phrase_atom',
                    2100037: 'phrase_atom',
                },
                'function': {
                    2000036: 'Conj',
                    2000037: 'Time',
                },
                'typ': {
                    2000036: 'CP',
                    2000037: 'PP',
                    2100036: 'CP',
                    2100037: 'PP',
                },
                'rela': {
                    2100036: 'NA',
                    2100037: 'NA',
                },
            },
            description='Ex 34:21; split up phrase',
        ),
        EditAction(
             edge_updates={
                'oslots': {
                    # phrases
                    695295: {
                        76399, 76400, 76401,
                    },
                    2000038: {
                        76402,
                    },
                    2000039: {
                        76403, 76404
                    },
                    2000040: {
                        76405,
                    },
                    2000041: {
                        76406, 76407, 76408
                    },
                    # phrase_atoms
                    950805: {
                        76399, 76400, 76401,
                    },
                    2100038: {
                        76402,
                    },
                    2100039: {
                        76403, 76404
                    },
                    2100040: {
                        76405,
                    },
                    2100041: {
                        76406, 76407, 76408
                    },
                },
                'head': {
                    76399: {695295},
                    76402: {2000038},
                    76403: {2000039},
                    76405: {2000040},
                    76406: {2000041},
                },
                'nhead': {
                    76400: {695295},
                    76404: {2000039},
                    76407: {2000041},
                },
             },
             feature_updates={
                'otype': {
                    2000038: 'phrase',
                    2000039: 'phrase',
                    2000040: 'phrase',
                    2000041: 'phrase',
                    2100038: 'phrase_atom',
                    2100039: 'phrase_atom',
                    2100040: 'phrase_atom',
                    2100041: 'phrase_atom',
                },
                'function': {
                    2000038: 'Conj',
                    2000039: 'Time',
                    2000040: 'Conj',
                    2000041: 'Time',
                },
                'typ': {
                    2000038: 'CP',
                    2000039: 'PP',
                    2000040: 'CP',
                    2000041: 'PP',
                    2100038: 'CP',
                    2100039: 'PP',
                    2100040: 'CP',
                    2100041: 'PP',
                },
                'rela': {
                    2100038: 'NA',
                    2100039: 'NA',
                    2100040: 'NA',
                    2100041: 'NA',
                },
             },
             description='Num 10:10; split up phrase'
        ),
    ],
    update_features={
        "function": {

            # Subj->Time
            # Gen 1:5; based on other wayehi-x clauses with
            # time phrases, this function appears to be mislabeled...
            # it's not that "evening was"...there is a dummy subject;
            # so: "*it* was evening"
            651617: "Time",

            # Subj-> Time
            # Gen 1:5; see note on 651617
            651620: "Time",

            # Adju->Time
            # Exod 12:18; mislabeled as Adju
            673370: "Time",

            # Adju->Time
            # Exod 12:18; mislabeled as Adju
            673373: "Time",

            # Adju->Time
            # Exod 12:18; mislabled as Adju
            673378: "Time",

            # Adju->Time
            # Deut 31:10; mislabled as Adju (חג)
            714507: "Time",

            # Time->Time
            # 1 Sam 25: 15; mislabled as Loca
            741508: "Time",

            # Time->Modi
            # 1 Sam 1:23; no reviewed sources took this as temporal
            731946: "Modi",

            # Time->PreC
            # Gen 9:29; phrase is pred. complement
            654059: "PreC",

            # Time->Adju
            # Josh 4:18; waters flowed 'as before'; not temporal
            716817: "Adju",

            # Time->Adju
            # 2 Kgs 13:5; 'as before'; not temporal
            769441: "Adju",

            # Time->Adju
            # 1 Sam 19:7; 'as before'; not temporal
            738950: "Adju",

            # Time->Adju
            # 1 Sam 21:6; 'as before'; not temporal
            739973: "Adju",

            # Time->Adju
            # Exod 30:10; This is frequentive not time location / duration
            679412: "Adju",

            # Time->Adju
            # Exod 30:10; This is frequentive not time location / duration
            679414: "Adju",

            # Time->Adju
            # Lev 16:34; This is frequentive not time location / duration
            688372: "Adju",

            # Time->Adju
            # Lev 23:16; how far one should count, not temporal
            690340: "Adju",

            # Time->Adju
            # Lev 26:18; No transs interprets as temporal
            691647: "Adju",

            # Time->Adju
            # Num 2:31; 'set out last', this is sequential but not with respect to timeline
            692794: "Adju",

            # Time->Cmpl
            # Josh 8:14; most take as locative; LXX omits; though Targum reads as ZMN, 'time'
            718111: "Cmpl",

            # Time->Adju
            # 1 Sam 18:10; comparative not temporal
            738537: "Adju",

            # Time->PreC
            # Num 28:14; the phrase belongs to prev phrase as part of genitive; 'חדש בחדו' is the
            # complete idiom; =burnt offfering of each month
            701713: "PreC",

            # Time->Freq
            # 1 Kgs 10:22; this frequentive explains some L+time cx
            757981: "Freq",

            # Time->Adju
            # 1 Kgs 7:24; NET does not take as temporal and LXX lacks this; probably not temporal
            756305: "Adju",

            # Time->Adju
            # Gen 48:7; this is a debatable example. NET does not read directly as temporal;
            # due to the difficulty, I'd like to exclude this from the sample set
            668696: "Adju",

            # Time->Freq
            # Ex 16:21; a repeated B-time is frequentive
            674998: "Freq",

            # Time->Freq
            # Lev 24:8
            690650: "Freq",

            # Time->Freq
            # 2 Sam 13:4
            747922: "Freq",

            # Time->Freq
            # Ex 30:7
            679386: "Freq",

            # Time->Freq
            # Ex 36:3
            681342: "Freq",

            # Time->Freq
            # Lev 6:5
            684164: "Freq",
        },
    },
    update_metadata={
        '': {
            "corpus": "BHSA-KinghamThesis",
            "description": "A modified version of the ETCBC's BHSA for my Cambridge PhD thesis",
            "version": "1.0",
            "editor": "Cody Kingham",
            "source": "Eep Talstra Centre for Bible and Computer",
            "source-url": "https://github.com/etcbc/bhsa",
            "encoders": "Constantijn Sikkel (QDF), Ulrik Petersen (MQL) and Dirk Roorda (TF)",
        },
        'omap@2021-KT': {
            'description': 'Mapping between nodes in BHSA 2021 version to BHSA Kingham Thesis version',
            'valueType': 'int',
        },
        'target': {
            **annotation_metadata,
            'description': 'Target value assigned by the thesis annotation project',
            'valueType': 'str',
        },
        'cl_type': {
            **annotation_metadata,
            'description': 'Clause type values assigned by the thesis annotation project',
            'valueType': 'str',
        },
        'aspect': {
            **annotation_metadata,
            'description': 'Clause aspect values assigned by the thesis annotation project',
            'valueType': 'str',
        },
        'tp_cluster': {
            **annotation_metadata,
            'description': 'Time phrase cluster values assigned by the thesis annotation project',
            'valueType': 'str',
        },
        'tense': {
            **annotation_metadata,
            'description': 'Clause tense values assigned by the thesis annotation project',
            'valueType': 'str',
        },
    },
)
