select
  *
from cv_endpoints.any_angina_prevalent_or_incident
full outer join (
  select sample_id, any_death
  from cv_endpoints.any_death
) any_death using (sample_id)

full outer join (
  select sample_id, cad_death
  from cv_endpoints.cad_death
) cad_death using (sample_id)

full outer join (
  select sample_id, cad_prevalent_or_incident
  from cv_endpoints.cad_prevalent_or_incident
) cad_prevalent_or_incident using (sample_id)

full outer join (
  select sample_id, cv_death
  from cv_endpoints.cv_death
) cv_death using (sample_id)

full outer join (
  select sample_id, mi_prevalent_or_incident
  from cv_endpoints.mi_prevalent_or_incident
) mi_prevalent_or_incident using (sample_id)

full outer join (
  select sample_id, pci_cabg_prevalent_or_incident
  from cv_endpoints.pci_cabg_prevalent_or_incident
) pci_cabg_prevalent_or_incident using (sample_id)

full outer join (
  select sample_id, unstable_angina_prevalent_or_incident
  from cv_endpoints.unstable_angina_prevalent_or_incident
) unstable_angina_prevalent_or_incident using (sample_id)

full outer join (
  select sample_id, hf as hf_prevalent_or_incident
  from legaultm_hcn4_paper.hf
) hf_prevalent_or_incident using (sample_id)
