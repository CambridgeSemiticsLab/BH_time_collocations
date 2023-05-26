"""Configure the build of BHSA KinghamThesis version."""

from corpus_building.utils import EditAction


# Corpus edits are logged here;
# NB: All node numbers are from ETCBC 2021 version
# New nodes will be assigned with the following scheme:
#     phrases = 2000000 + N
#     phrase_atoms = 2001000 + N
#     clauses = 3000000 + N
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
                }
            },
            deletions={661279},
            description="Gen 30:16; BLJLH H>W, merge phrases",
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    662567: {17889, 17890, 17891, 17892},
                },
            },
            deletions={662568},
            description="Gen 32:23;  BLJLH H>W, merge phrases",
        ),
        EditAction(
            edge_updates={
                'oslots': {
                    435277: {
                        37530, 37531, 37532, 37533, 37534,
                        37535, 37536, 37537, 37540, 37541,
                        37542, 37543, 37544
                    },
                },
            },
            deletions={435279},
            description="Ex 16:8; merge elliptical clause into prev",
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
                },
                'mother': {
                    929211: {929210},  # phrase_atom
                    523704: {523702},  # clause_atom
                },
            },
            deletions={435374, 523703},
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
                    690894: {
                        67313, 67314, 67315, 67316, 67317
                    },
                    2000005: {
                        67318, 67319, 67320, 67321, 67322, 67323
                    },
                },
                'head': {67318: {2000005}},
                'nhead': {67320: {2000005}},
            },
            feature_updates={
                'otype': {2000005: 'phrase'},
                'function': {2000005: 'Time'},
                'typ': {2000005: 'PP'},
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
                    2001001: {82442},
                    2001002: {
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
                    2001001: 'phrase_atom',
                    2001002: 'phrase_atom',
                },
                'function': {
                    2000012: 'Conj',
                    2000013: 'Time',
                },
                'typ': {
                    2000012: 'CP',
                    2000013: 'PP',
                    2001001: 'CP',
                    2001002: 'PP',
                },
                'rela': {
                    2001001: 'NA',
                    2001002: 'NA',
                },
            },
            description='Num 9:19; split up phrase',
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
    },
)
