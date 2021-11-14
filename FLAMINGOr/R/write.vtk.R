#' write.vtk
#'
#' Convert the FLAMINGO predicted 3D genome structure into .vtk format for visualization
#' @param points 3D coordiantes predicted by FLAMINGO in x,y,z format.
#' @param lookup_table The annotation for each point, could be labels or scores, i.e. TAD annotations.
#' @param name output file name
#' @param opt_path output file path
#' @keywords FLAMINGO
#' @return Write out a .vtk file for visualization using Paraview
#' @export
write.vtk <- function(points,lookup_table,name,opt_path){
  write.table("# vtk DataFile Version 1.0",opt_path,col.names = F,row.names = F,sep = "\n",quote=F)
  # cat("# vtk DataFile Version 1.0\n")
  write.table(name,opt_path,col.names = F,row.names = F,sep = "\n",quote=F,append = T)
  write.table("ASCII",opt_path,col.names = F,row.names = F,sep = "\n",quote=F,append = T)
  write.table("\nDATASET POLYDATA",opt_path,col.names = F,row.names = F,sep = "\n",quote=F,append = T)
  n_points <- dim(points)[1]
  write.table(paste("POINTS",n_points,"float"),opt_path,col.names = F,row.names = F,sep = "\n",quote=F,append = T)
  write.table(points,opt_path,col.names = F,row.names = F,quote=F,append = T)
  write.table("\n",opt_path,col.names = F,row.names = F,sep = "\n",quote=F,append = T)
  write.table(paste("LINES 1",n_points+1,n_points),opt_path,col.names = F,row.names = F,sep = "\n",quote=F,append = T)
  lines <- 0:(n_points-1)
  write.table(lines,opt_path,col.names = F,row.names = F,sep = "\n",quote=F,append = T)
  write.table("\n",opt_path,col.names = F,row.names = F,sep = "\n",quote=F,append = T)
  write.table(paste("POINT_DATA",n_points),opt_path,col.names = F,row.names = F,sep = "\n",quote=F,append = T)
  write.table(paste("SCALARS volume float\n","LOOKUP_TABLE default"),opt_path,col.names = F,row.names = F,sep = "\n",quote=F,append = T)
  write.table(lookup_table,opt_path,col.names = F,row.names = F,sep = "\n",quote=F,append = T)
}
