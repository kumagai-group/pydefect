from typing import List

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.calc_summary import CalcSummary, SingleCalcSummary
from pydefect.analyzer.make_energy_summary import num_atom_differences
from pydefect.defaults import defaults
from pydefect.input_maker.defect_entry import DefectEntry


def make_calc_summary(defect_entry_list: List[DefectEntry],
                      calc_results_list: List[CalcResults],
                      p_calc_results: CalcResults) -> CalcSummary:
    summaries = {}
    for entry, result in zip(defect_entry_list, calc_results_list):
        atom_io = num_atom_differences(
            result.structure, p_calc_results.structure, str_key=True)
        relative_energy = result.energy - p_calc_results.energy
        is_energy_strange = abs(relative_energy) > defaults.abs_strange_energy
        summaries[entry.full_name] = SingleCalcSummary(entry.charge,
                                                       atom_io,
                                                       result.electronic_conv,
                                                       result.ionic_conv,
                                                       is_energy_strange)
    return CalcSummary(single_summaries=summaries)
