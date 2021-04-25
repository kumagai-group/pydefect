from typing import List

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.calc_summary import CalcSummary, SingleCalcSummary
from pydefect.analyzer.defect_structure_info import DefectStructureInfo
from pydefect.analyzer.make_defect_energy_info import num_atom_differences
from pydefect.defaults import defaults
from pydefect.input_maker.defect_entry import DefectEntry


def make_calc_summary(defect_entry_list: List[DefectEntry],
                      calc_results_list: List[CalcResults],
                      structure_list: List[DefectStructureInfo],
                      p_calc_results: CalcResults) -> CalcSummary:
    summaries = {}
    for entry, result, str_info in zip(defect_entry_list, calc_results_list,
                                       structure_list):
        atom_io = num_atom_differences(result.structure,
                                       p_calc_results.structure)
        relative_energy = result.energy - p_calc_results.energy
        is_energy_strange = abs(relative_energy) > defaults.abs_strange_energy
        summaries[entry.full_name] = SingleCalcSummary(
            charge=entry.charge,
            atom_io=atom_io,
            electronic_conv=result.electronic_conv,
            ionic_conv=result.ionic_conv,
            is_energy_strange=is_energy_strange,
            same_config_from_init=str_info.same_config_from_init,
            defect_type=str_info.defect_type,
            symm_relation=str_info.symm_relation)
    return CalcSummary(single_summaries=summaries)
