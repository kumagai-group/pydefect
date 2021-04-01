# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.


user_incar_settings = {"LPARD": True,
                       "LSEPB": True,
                       "KPAR": 1,
                       "IBAND": args.band_indices}

if args.kpoint_indices:
    user_incar_settings["KPUSE"] = args.kpoint_indices

vasp_set = ViseInputSet.from_prev_calc(
    dirname=args.read_dir,
    parse_calc_results=False,
    parse_incar=True,
    sort_structure=False,
    standardize_structure=False,
    files_to_transfer={"WAVECAR": "L"},
    user_incar_settings=user_incar_settings,
    contcar_filename=args.contcar)

vasp_set.write_input(args.write_dir)