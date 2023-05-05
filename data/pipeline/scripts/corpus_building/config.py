"""Configure the build of BHSA Kingham Thesis version."""

from corpus_building.utils import EditAction


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
        # NB: below is just an example
        # EditAction(
        #     feature_updates={
        #         'function': {974852: 'LocaTime'},
        #         'prep_type': {974852: 'B_simul'},
        #     },
        #     edge_updates={
        #         'oslots': {
        #             484384: {283995, 283996, 283997, 283998,
        #                      283999, 284000, 284001},
        #             484385: {284002, 284003, 284004, 284005},
        #         },
        #     },
        # ),
    ],
    update_features={
        "function": {
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
            "corpus": "BHSA-Kingham-thesis",
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
        'prep_type': {'description': 'test123', 'valueType': 'str'},
    },
)
