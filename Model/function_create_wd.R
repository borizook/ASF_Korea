create_wd = function(folder_scen){
  
  if (folder_scen %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_scen))}
  
  scen_folder_wd = file.path(getwd(), folder_scen)
  
  setwd(scen_folder_wd)
  
  folder_init_data = "init_data"
  if (folder_init_data %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_init_data))}
  
  
  folder_pars_used = "Pars_used"
  if (folder_pars_used %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_pars_used))}
  
  
  folder_init_data = "init_data"
  if (folder_init_data %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_init_data))}
  
  folder_artificial_dataset = "Artificial_dataset"
  if (folder_artificial_dataset %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_artificial_dataset))}
  
  folder_diag = "Diagnostic"
  if (folder_diag %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_diag))}
  
  
  folder_traj = "Trajectories"
  if (folder_traj %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_traj))}
  
  
  folder_traj_sp = "Trajectories/Spatial"
  if (folder_traj_sp %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_traj_sp))}
  
  
  
  folder_traj_tem = "Trajectories/Temporal"
  if (folder_traj_tem %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_traj_tem))}
  
  
  folder_post = "Posterior distribution"
  if (folder_post %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_post))}
  
  setwd(file.path(getwd(), folder_post))
  
  folder_overlap = "overlap"
  
  if (folder_overlap %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_overlap))}
  
  
  folder_dist_plot = "dist_plot"
  
  if (folder_dist_plot %in% dir(getwd()) == FALSE) {
    dir.create(file.path(getwd(), folder_dist_plot))}
  
  scen_folder_wd <<- scen_folder_wd
  folder_artificial_dataset <<- folder_artificial_dataset
  folder_init_data <<- folder_init_data
  folder_post<<-folder_post
  folder_traj_tem<<-folder_traj_tem
  folder_traj_sp<<- folder_traj_sp
  folder_traj<<- folder_traj
  folder_diag<<- folder_diag
  folder_artificial_dataset<<- folder_artificial_dataset
  folder_init_data<<- folder_init_data
  folder_pars_used<<- folder_pars_used
  setwd(scen_folder_wd)
  
  
}
