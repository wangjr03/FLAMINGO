#' flamingo.assemble_structure
#'
#' Assemble the domain-level structure into the genome-wide structure
#' @param flamingo_backbone_prediction_obj A flamingo_prediction object containing the 3D structure of the backbone
#' @param flamingo_final_res_data_obj A flamingo object containing the IF and PD matrix under the fine resolution, e.g. 5kb.
#' @param list_of_flamingo_domain_prediction_obj A list of flamingo_prediction object containing the intra-domain 3D structures.
#' @param max_iter maximum number of iteration, default is 500.
#' @keywords FLAMINGO
#' @return A data.frame containing the FLAMINGO predicted 3D structure.
#' @export

flamingo.assemble_structure <- function(flamingo_backbone_prediction_obj,flamingo_final_res_data_obj,list_of_flamingo_domain_prediction_obj,max_iter=500){
  n = dim(flamingo_backbone_prediction_obj@coordinates)[1]
  backbone_id <- flamingo_backbone_prediction_obj@id
  backbone <- as.matrix(flamingo_backbone_prediction_obj@coordinates)
  n_point_per_domain = list_of_flamingo_domain_prediction_obj[[1]]@input_n
  all_points <- list()
  counter_idx <- 1
  rm_id_list <- list()
  val_backbond_id <- c()
  scaler <- as.matrix(dist(backbone))
  scaler <- mean(scaler[which(row(scaler)+1 == col(scaler))])
  print('Reading substructures...')
  for(counter in backbone_id){
    if(!counter %in% names(list_of_flamingo_domain_prediction_obj)){
      next
    }
    tmp_domain_res = list_of_flamingo_domain_prediction_obj[[as.character(counter)]]
    tmp_id <- tmp_domain_res@id
    p_t <- as.matrix(tmp_domain_res@coordinates)[tmp_id,]
    tmp_center <- backbone[which(backbone_id==counter),]
    rm_id <- setdiff(1:n,tmp_id)
    rm_id_list[[counter_idx]] <- rm_id
    
    if(dim(p_t)[1] == 0){
      next
    }
    old_center <- apply(p_t,2,mean)
    radius <- max(apply(p_t,1,function(x){norm(x-old_center,'2')}))
    new_loc <- t(apply(p_t,1,function(x){
      
      (x-old_center)/radius*scaler+tmp_center
      
    }))
    all_points[[counter_idx]] <- as.data.frame(new_loc)
    counter_idx = counter_idx+1
    val_backbond_id <- c(val_backbond_id,counter)
  }
  frag_id <- val_backbond_id
  names(all_points) <- frag_id
  total_struc <- backbone[match(val_backbond_id,backbone_id),]
  #prepare input data
  #calculate point id
  id_list <- list()
  for(i in 1:length(frag_id)){
    tmp_id <- 1:n_point_per_domain
    if(length(rm_id_list[[i]])>0){
      tmp_id <- tmp_id[-rm_id_list[[i]]]
    }
    
    tmp_id <- (frag_id[i]-1)*n_point_per_domain+tmp_id
    id_list[[i]] <- tmp_id[1:dim(all_points[[i]])[1]]
  }
  
  direction <- (total_struc[2,]-total_struc[1,])
  tail_direction <- tail(all_points[[1]],n=1)-total_struc[1,]
  direction <- direction/norm(direction,'2')
  tail_direction <- tail_direction/norm(tail_direction,'2')
  direction <- as.vector(direction)
  tail_direction <- as.vector(as.matrix(tail_direction))
  r_mat <- rotation_matrix(matrix(tail_direction,ncol=1),as.matrix(direction,nrow=1))
  all_points[[1]] <- rotate(all_points[[1]],total_struc[1,],r_mat)
  pd <- flamingo_final_res_data_obj@PD
  error <- 10
  d = get_dist_vec(all_points,id_list,pd)
  t = 1e-4
  #all_points_bp <- all_points
  g_p <- matrix(0,length(all_points)-1,3)
  y_s_p <-  matrix(0,length(all_points)-1,3)
  y_e_p <-  matrix(0,length(all_points)-1,3)
  error_p = error_change <- 1
  N <- dim(do.call(rbind,all_points))[1]
  p_p <- matrix(0,N,3)
  iter <- 1
  while(error>1e-4&error_change>1e-5){
    y_s <- get_start_point(all_points)
    y_e <- get_end_point(all_points)
    d_tild <- evaluate_dist(y_s,y_e)
    g <- 4*(d_tild-d)*(y_s-y_e)
    
    y = g-g_p
    s = y_s-y_s_p
    t_k = (sum(diag(t(s)%*%s)))/(sum(diag(t(s)%*%y)))
    
    y_s_prim <- y_s-t_k*g
    g_p <- g
    y_s_p <- y_s
    #update structure
    for(i in 2:length(all_points)){
      tmp_points <- all_points[[i]]
      tmp_start <- tmp_points[1,]
      tmp_center <- total_struc[i,]
      old_direction <- as.vector(as.matrix(tmp_start-tmp_center))+1e-3
      new_direction <- y_s_prim[i-1,] - tmp_center
      r_mat <- rotation_matrix(old_direction,new_direction)
      all_points[[i]] <- rotate(all_points[[i]],total_struc[i,],r_mat)
    }
    #update structure
    if(iter == 1){
      y_s <- get_start_point(all_points)
      y_e <- get_end_point(all_points)
      d_tild <- evaluate_dist(y_s,y_e)
      g <- 4*(d_tild-d)*(y_s-y_e)
      y = g-g_p
      s = y_e-y_e_p
      t_k = (sum(diag(t(s)%*%s)))/(sum(diag(t(s)%*%y)))
      y_e_prim <- y_e-t_k*g
      y_e_p <- y_e
      for(i in 1:(length(all_points)-1)){
        tmp_points <- all_points[[i]]
        tmp_end <- tail(tmp_points,n=1)
        tmp_center <- total_struc[i,]
        old_direction <- as.vector(as.matrix(tmp_end-tmp_center))+1e-3
        new_direction <- y_e_prim[i,] - tmp_center
        r_mat <- rotation_matrix(old_direction,new_direction)
        all_points[[i]] <- rotate(all_points[[i]],total_struc[i,],r_mat)
      }
    }
    iter <- iter + 1
    cur_p <- do.call(rbind,all_points)
    all_points_err <- norm(cur_p-p_p,'F')/norm(cur_p,'F')
    p_p <- cur_p
    error <- norm(d_tild-d,'2')/norm(d,'2')
    error_change <- abs(error-error_p)
    error_p <- error
    print(error)
    if(all_points_err<0.0001){
      break
    }
    if(iter > max_iter){
      break
    }
    
  }
  all_points <- smt(do.call(rbind,all_points))
  res = data.frame(frag_id = unlist(id_list),x=all_points[,1],y=all_points[,2],z=all_points[,3])
  return(res)
}
