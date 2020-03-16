get.cell.idx.within.dist.range <- function(point.ref, points, r.i, r.f, index = TRUE, dim = 3){
  #point.ref <- as.numeric(data.tconv.highest[1,])
  #points <- as.matrix(data.treg.Nur)
  num.row <- nrow( points)
  
  point.ref = t(matrix(rep(point.ref, num.row ), ncol =num.row ))
  if(dim == 3){
    dist.euc <- points - point.ref
    dist.euc <- sqrt(rowSums(dist.euc^2))
    # r.i <- 2
    # r.f <- 3
    idx <- dist.euc >= r.i & dist.euc < r.f
    if(index){
      return(idx)
    }else{
      if(is.null(nrow(points[idx,]))){
        return(cbind(t(points[idx,]),dist = dist.euc[idx], idx = (1:num.row)[idx]))
      }else{
        return(cbind(points[idx,],dist = dist.euc[idx], idx = (1:num.row)[idx]))
      }
    }
  }else if(dim == 2){
    dist.euc <- points[,1:2] - point.ref[,1:2]
    dist.euc <- sqrt(rowSums(dist.euc^2))
    # r.i <- 2
    # r.f <- 3
    idx <- dist.euc >= r.i & dist.euc < r.f
    if(index){
      return(idx)
    }else{
      if(is.null(nrow(points[idx,]))){
        return(cbind(t(points[idx,]),dist = dist.euc[idx], idx = (1:num.row)[idx]))
      }else{
        return(cbind(points[idx,],dist = dist.euc[idx], idx = (1:num.row)[idx]))
      }
    }
  }
}


extract.cells.increasing.shells <- function(cells.list, refs.list, shells){
  r.i <- shells
  idx.list <- names(cells.list)
  print(idx.list)
  cells.extract <- list()
  
  for(ln in idx.list){ ## for each LN section
    print(ln)
    cells.extract[[ln]] <- data.frame()
   # print(names(cells.extract))
    
    for(r in r.i){
      print(r)
      for(i in 1:nrow(refs.list[[ln]]) ){ ## for ref T cells and total tregs
        temp <- get.cell.idx.within.dist.range(
          point.ref = as.numeric(refs.list[[ln]][i,1:3]),
          points = as.matrix(cells.list[[ln]][,1:3]),
          r.i = r,
          r.f = r+1,
          index = FALSE)
       # print(nrow(temp))
        if(nrow(temp)!=0){
          num.col <- ncol(cells.list[[ln]])
          #print(num.col)
          #print(dim(temp))
          if(num.col >3){
            cells.extract[[ln]] <- rbind(cells.extract[[ln]],                              
                                         cbind(temp, 
                                               center.idx = i,
                                               r.i = r,
                                               r.f = r+1,
                                               cells.list[[ln]][temp[,5],4:num.col])
            )
          }else{
            cells.extract[[ln]] <- rbind(cells.extract[[ln]],                              
                                         cbind(temp, 
                                               center.idx = i,
                                               r.i = r,
                                               r.f = r+1
                                               )
            )
          }
                                         
        
      
        }
      }
    }
  }
  return(cells.extract)
}


extract.cells.within.dist <- function(cells.list, refs.list, distance){
  idx.list <- names(cells.list)
  cells.extract <- list()
  
  for(ln in idx.list){ ## for each LN section
    
    cells.extract[[ln]] <- NULL
    
    
    #for(r in r.i){
      for(i in 1:nrow(refs.list[[ln]]) ){ ## for ref T cells and total tregs
        print(i)
        temp <- get.cell.idx.within.dist.range(
          point.ref = as.numeric(refs.list[[ln]][i,1:3]),
          points = as.matrix(cells.list[[ln]][,1:3]),
          r.i = 0,
          r.f = distance+1,
          index = FALSE)
        if(nrow(temp)!=0){
          num.col <- ncol(cells.list[[ln]])
          #print(num.col)
          #print(dim(temp))
          if(num.col >3){
            cells.extract[[ln]] <- rbind(cells.extract[[ln]],                              
                                         cbind(temp, 
                                               center.idx = i,
                                               r.i = 0,
                                               r.f = distance+1,
                                               cells.list[[ln]][temp[,5],4:num.col])
            )
          }else{
            cells.extract[[ln]] <- rbind(cells.extract[[ln]],                              
                                         cbind(temp, 
                                               center.idx = i,
                                               r.i = 0,
                                               r.f = distance+1
                                         )
            )
          }
          
          
          
        }
      }
    #}
  }
  return(cells.extract)
}


#count the number of cells per shell

count.number.per.shell <- function(extracted.list, shells){
  num.per.shell <- list()
  num.cumul <- list()
  idx.list <- names(extracted.list)
  #print(idx.list)
  for(ln in idx.list){
    print(ln)
    num.cumul[[ln]] <- data.frame()
    num.per.shell[[ln]] <- data.frame()

    dist.temp <- as.character(shells)
    
    if(nrow(extracted.list[[ln]]) != 0){
      extracted.list[[ln]] <- data.frame(extracted.list[[ln]])
      ##total Tregs
      length.center.idx <- as.numeric(names(table(extracted.list[[ln]][,"center.idx"])))
      for(center.idx in length.center.idx){
        temp.count <- NULL
        temp.count.cumul <- NULL
        #temp.count <- append(temp.count,center.idx)
        temp.table <- table(extracted.list[[ln]]$r.i[extracted.list[[ln]]$center.idx ==center.idx])
        for(i in dist.temp){
          if(is.na(temp.table[i])){
            temp.count <- append(temp.count,0)
          }else{
            temp.count <- append(temp.count,temp.table[i] )
          }
        }
        #print("a")
        temp.count.cumul <- append(temp.count.cumul, temp.count[1])
        for(i in 2:length(shells)){
          temp.count.cumul <-append(temp.count.cumul, temp.count[i] +temp.count.cumul[i-1])
          
        }
        
        temp.count
        temp.count.cumul
        
        num.cumul[[ln]] <- rbind(num.cumul[[ln]], append(center.idx,temp.count.cumul))
        num.per.shell[[ln]] <- rbind(num.per.shell[[ln]], append(center.idx,temp.count))
        
      }
      #print("b")
      names( num.cumul[[ln]]) <- append("center.idx", shells + shells[2]-shells[1])
      names( num.per.shell[[ln]] ) <- append("center.idx", shells + shells[2]-shells[1])
    }
  
  }
  return(list(num.cumul = num.cumul, num.per.shell = num.per.shell))
}

count.number.within.dist <- function(extracted.list){
  num.per.shell <- list()
  num.cumul <- list()
  idx.list <- names(extracted.list)
  #print(idx.list)
  for(ln in idx.list){
    print(ln)
    num.cumul[[ln]] <- data.frame()
   # num.per.shell[[ln]] <- data.frame()
    
    #dist.temp <- as.character(shells)
    
    extracted.list[[ln]] <- data.frame(extracted.list[[ln]])
    
    num.rows <- nrow(extracted.list[[ln]])
    ##total Tregs
    length.center.idx <- as.numeric(names(table(extracted.list[[ln]][,"center.idx"])))
    for(center.idx in length.center.idx){
      temp.count <- NULL
      #temp.count <- append(temp.count,center.idx)
      temp.count <- num.rows - table(extracted.list[[ln]]$center.idx ==center.idx)["FALSE"]
      
      #print("a")
      #temp.count
      num.cumul[[ln]] <- rbind(num.cumul[[ln]], append(center.idx,temp.count))
    }
    #print("b")
    names( num.cumul[[ln]]) <- append("center.idx", "num.cells")
  }
  return(num.cumul)
}


extract.count.cells.within.dist <- function(cells.list, refs.list, distance){
  idx.list <- names(cells.list)
  num.cumul <- list()
  for(ln in idx.list){ ## for each LN section
    print(ln)
    num.cumul[[ln]] <- data.frame()
    for(i in 1:nrow(refs.list[[ln]]) ){ ## for ref T cells and total tregs
      if(i%%1000 ==0){
        print(i)
      }
      temp <- get.cell.idx.within.dist.range(
        point.ref = as.numeric(refs.list[[ln]][i,1:3]),
        points = as.matrix(cells.list[[ln]][,1:3]),
        r.i = 0,
        r.f = distance+1,
        index = FALSE)
      temp.count <- NULL
      temp.count <- nrow(temp)
      #print(temp)
      if(temp.count!=0){
        #print("a")
        #temp.count
        #temp.row <-  rbind(num.cumul[[ln]], data.frame(i,temp.count,paste0(temp[,5], collapse = ","),paste0(temp[,4], collapse = ",") ))
        #print(nrow(temp.row ))
        num.cumul[[ln]] <- rbind(num.cumul[[ln]], data.frame(center.idx = i,
                                                             num.cells = temp.count,
                                                             cell.idx = paste0(temp[,5], collapse = ","),
                                                             cell.dist = paste0(temp[,4], collapse = ","),
                                                             stringsAsFactors = F))
        #print(num.col)
        #print(dim(temp))
      }else{
        num.cumul[[ln]] <- rbind(num.cumul[[ln]], data.frame(center.idx = i,
                                                             num.cells = temp.count,
                                                             cell.idx = "",
                                                             cell.dist = "",
                                                             stringsAsFactors = F ))
      }
    }
  }
  return(num.cumul)
}


extract.cells.increasing.shells <- function(cells.list, refs.list, shells){
  r.i <- shells
  idx.list <- names(cells.list)
  print(idx.list)
  cells.extract <- list()
  
  for(ln in idx.list){ ## for each LN section
    print(ln)
    cells.extract[[ln]] <- data.frame()
    # print(names(cells.extract))
    
    for(r in r.i){
      print(r)
      for(i in 1:nrow(refs.list[[ln]]) ){ ## for ref T cells and total tregs
        temp <- get.cell.idx.within.dist.range(
          point.ref = as.numeric(refs.list[[ln]][i,1:3]),
          points = as.matrix(cells.list[[ln]][,1:3]),
          r.i = r,
          r.f = r+1,
          index = FALSE)
        # print(nrow(temp))
        if(nrow(temp)!=0){
          num.col <- ncol(cells.list[[ln]])
          #print(num.col)
          #print(dim(temp))
          if(num.col >3){
            cells.extract[[ln]] <- rbind(cells.extract[[ln]],                              
                                         cbind(temp, 
                                               center.idx = i,
                                               r.i = r,
                                               r.f = r+1,
                                               cells.list[[ln]][temp[,5],4:num.col])
            )
          }else{
            cells.extract[[ln]] <- rbind(cells.extract[[ln]],                              
                                         cbind(temp, 
                                               center.idx = i,
                                               r.i = r,
                                               r.f = r+1
                                         )
            )
          }
          
          
          
        }
      }
    }
  }
  return(cells.extract)
}


extract.count.cells.shell <- function(cells.list, refs.list, shells){
  r.i <- shells
  idx.list <- names(cells.list)
  num.cumul <- list()
  for(ln in idx.list){ ## for each LN section
    print(ln)
    num.cumul[[ln]] <- data.frame()
    for(r in r.i){
      print(r)
      for(i in 1:nrow(refs.list[[ln]]) ){ ## for ref T cells and total tregs
        if(i%%1000 ==0){
          print(i)
        }
        temp <- get.cell.idx.within.dist.range(
          point.ref = as.numeric(refs.list[[ln]][i,1:3]),
          points = as.matrix(cells.list[[ln]][,1:3]),
          r.i = r,
          r.f = r+1,
          index = FALSE)
        temp.count <- NULL
        temp.count <- nrow(temp)
        #print(temp)
        if(temp.count!=0){
          #print("a")
          #temp.count
          #temp.row <-  rbind(num.cumul[[ln]], data.frame(i,temp.count,paste0(temp[,5], collapse = ","),paste0(temp[,4], collapse = ",") ))
          #print(nrow(temp.row ))
          num.cumul[[ln]] <- rbind(num.cumul[[ln]], data.frame(center.idx = i,
                                                               r.i = r,
                                                               num.cells = temp.count,
                                                               cell.idx = paste0(temp[,5], collapse = ","),
                                                               cell.dist = paste0(temp[,4], collapse = ","),
                                                               stringsAsFactors = F))
          #print(num.col)
          #print(dim(temp))
        }else{
          num.cumul[[ln]] <- rbind(num.cumul[[ln]], data.frame(center.idx = i,
                                                               r.i = r,
                                                               num.cells = temp.count,
                                                               cell.idx = "",
                                                               cell.dist = "",
                                                               stringsAsFactors = F ))
        }
      }
    }
  }
  return(num.cumul)
}


count.number.per.shell.1 <- function(tregs.extract.shell){
  tregs.extract.num.1um <- list()
  
  for(idx.list in names(tregs.extract.shell)){
    center.idx <- unique(tregs.extract.shell[[idx.list]]$center.idx)
    tregs.extract.num.1um[[idx.list]] <- data.frame(center.idx = center.idx)
    r.f <- unique(tregs.extract.shell[[idx.list]]$r.i) + 1
    temp.mat <- NULL
    for(idx.cell in center.idx){
      temp.mat <- rbind(temp.mat, tregs.extract.shell[[idx.list]]$num.cells[tregs.extract.shell[[idx.list]]$center.idx==idx.cell])
    }
    tregs.extract.num.1um[[idx.list]] <- cbind(tregs.extract.num.1um[[idx.list]], temp.mat) 
  }
  temp.mat <- NULL
  for(idx.list in names(tregs.extract.shell)){
    temp.mat<- rbind(temp.mat,tregs.extract.num.1um[[idx.list]])
  }
  return(tregs.extract.num.1um )
  
}



sum.marker.intensity.per.shell <- function(extracted.list ,shells){
  sum.per.shell <- list()
  idx.list <- names(extracted.list)
  #print(idx.list)
  for(ln in idx.list){
    print(ln)
    sum.per.shell[[ln]] <- list()
    dist.temp <- as.character(shells)
    extracted.list[[ln]] <- data.frame(extracted.list[[ln]])
    num.col = ncol(extracted.list[[ln]])
    markers <- names(  extracted.list[[ln]] )[9:num.col]
    for(idx.marker in markers){
      sum.per.shell[[ln]][[idx.marker]] <- data.frame()
    }
    ##total Tregs
    length.center.idx <- length(table(extracted.list[[ln]][,"center.idx"]))
    for(center.idx in 1:length.center.idx){
      temp.sum <- list()
      for(idx.marker in markers){
        temp.sum[[idx.marker]] <- NULL
      }
      #temp.count <- append(temp.count,center.idx)
      #temp.table <- table(extracted.list[[ln]]$r.i[extracted.list[[ln]]$center.idx ==center.idx])
      temp.frame <- extracted.list[[ln]][extracted.list[[ln]]$center.idx ==center.idx ,c('r.i', markers)]
      for(i in dist.temp){
        # if(any(temp.frame$r.i==i)){
        #   for(idx.marker in markers){
        #     temp.sum[[idx.marker]] <- append(temp.sum[[idx.marker]], sum(temp.frame[temp.frame$r.i==i,idx.marker]))
        #   }
        # }else{
        #   for(idx.marker in markers){
        #     temp.sum[[idx.marker]] <- append(temp.sum[[idx.marker]], NA)
        #   }
        # }
        
        for(idx.marker in markers){
          temp.sum[[idx.marker]] <- append(temp.sum[[idx.marker]], sum(temp.frame[temp.frame$r.i==i,idx.marker]))
        }
      }
      #print("a")
   
      
      for(idx.marker in markers){
        sum.per.shell[[ln]][[idx.marker]] <-rbind(sum.per.shell[[ln]][[idx.marker]],append(center.idx,temp.sum[[idx.marker]]))
      }
      
    }
    
    for(idx.marker in markers){
      names(sum.per.shell[[ln]][[idx.marker]]) <-  append("center.idx", shells + shells[2]-shells[1])
    }
  }
  return(sum.per.shell)
}






if(1){#local machine
  Sys.setenv("PKG_CXXFLAGS"="-std=c++11 -I/Users/kyemyungpark/Dropbox/Codes/project_tcell_activation/analysis/libigl/include  -I/opt/local/include -I/opt/local/include/eigen3",
             "PKG_LIBS"="-lm -lCGAL -lCGAL_Core -frounding-math -lmpfr -lgmp ")
  sourceCpp("~/Dropbox/Codes/project_tcell_activation/analysis/mesh_boolean_R.cpp", rebuild = F, cacheDir = "~/Dropbox/Codes/project_tcell_activation/analysis/sharedlib/" )
  Sys.setenv("PKG_CXXFLAGS"="-std=c++11 -I/Users/kyemyungpark/Dropbox/Codes/project_tcell_activation/analysis/libigl/include -I/usr/local -I/opt/local/include -I/opt/local/include/eigen3", "PKG_LIBS"="-lm -lCGAL_Core -frounding-math -lmpfr -lgmp")
  sourceCpp("~/Dropbox/Codes/project_tcell_activation/analysis/mesh_volume_R.cpp", rebuild = , cacheDir = "~/Dropbox/Codes/project_tcell_activation/analysis/sharedlib/" )
  
}






if(0){
  #in cluster locus
  library(Rcpp)
  Sys.setenv("PKG_CXXFLAGS"="-std=c++11",
             "PKG_LIBS"="-lm -lCGAL -lCGAL_Core -frounding-math -lmpfr -lgmp -lboost_thread")
  sourceCpp("mesh_boolean_R.cpp", rebuild = T, verbose= T, cacheDir = "sharedlib/" )
  # -I/nethome/parkk6/.conda/envs/volumetry/include
  Sys.setenv("PKG_CXXFLAGS"="-std=c++11",
             "PKG_LIBS"="-lm -lCGAL -lCGAL_Core -frounding-math -lmpfr -lgmp -lboost_thread")
  sourceCpp("mesh_volume_R.cpp", rebuild = T, cacheDir = "sharedlib/" )
  
}


wrap_mesh_boolean <- function(mesh1, mesh2, op_type = "intersect" ){
  temp.output <- mesh_boolean(VAi = as.matrix(t(mesh1$vb[1:3,])),
                              FAi = as.matrix(t(mesh1$it)-1),
                              VBi = as.matrix(t(mesh2$vb[1:3,])),
                              FBi = as.matrix(t(mesh2$it)-1),
                              op_type = op_type)
  mesh_output <- list(
    vb = rbind(t(temp.output$VC), 1),
    it = t(temp.output$FC+1),
    primitivetype = "triangle",
    material = NULL,
    normals = NULL,
    texcoords = NULL
  )
  class(mesh_output) <- c("mesh3d",  "shape3d") 
  
  return(mesh_output)
}





#wrapper function of generating sphere

const.sphere.mesh <- function(r, center, subdivision = 4){
  #library(Rvcg)
  sphere.obj <-vcgSphere(subdivision = subdivision)
  #sphere.obj <- t(sphere.obj$vb[1:3,])
  sphere.obj$vb[1:3,] <-sphere.obj$vb[1:3,]*r # scale
  sphere.obj$vb[1:3,] <- sphere.obj$vb[1:3,] + matrix(rep(as.numeric(center),ncol(sphere.obj$vb)),nrow = 3)
  return(sphere.obj)
}


const.ellipse.mesh <- function(r, h, center, subdivision = 4){
  #library(Rvcg)
  sphere.obj <-vcgSphere(subdivision = subdivision)
  #sphere.obj <- t(sphere.obj$vb[1:3,])
  sphere.obj$vb[1:2,] <-sphere.obj$vb[1:2,]*r # scale
  sphere.obj$vb[3,] <-sphere.obj$vb[3,]*h # scale
  
  sphere.obj$vb[1:3,] <- sphere.obj$vb[1:3,] + matrix(rep(as.numeric(center),ncol(sphere.obj$vb)),nrow = 3)
  return(sphere.obj)
}



as.ashape3d <- function(mesh, alpha = 200, pert = F){
  return(ashape3d(t(mesh$vb[1:3,]), alpha = alpha, pert = pert))
}


const.cylinder.mesh <- function(r, height= 20, center){
  height = 20
  center = c(0,0,0)
  r= 1 
  center.top = center
  center.top[3] = center.top[3] + height
  center.bottom = center
  center.bottom[3] = center.bottom[3] - height
  # center = rbind(center.bottom, # start
  #                center.top ), 
 
  c= cylinder3d(center = rbind(center.bottom, # start
                               center.top ), # end
               radius = r,
               sides=100
  )
  
  
#  c <- subdivision3d(c, depth = subdivision)
  #plot3d(c)
  return(c)
  
}






wrap_mesh_volume <- function(mesh){
  temp.output <- mesh_volume(Vi = as.matrix(t(mesh$vb[1:3,])),
                             Fi = as.matrix(t(mesh$it)-1)
                            )
  return( temp.output)
}








##Defining surfaces
def.surface <- function(extracted.list, refs.list, alpha = 200){
  surface.list <- list()
  idx.list <- names(refs.list)
  for(ln in idx.list){
    print(ln)
    surface.list[[ln]] <-list()
    extracted.list[[ln]] <- data.frame(extracted.list[[ln]])
    for(i in 1:nrow(refs.list[[ln]])){
      temp.ashape <- ashape3d(as.matrix(extracted.list[[ln]][extracted.list[[ln]]$center.idx == i ,c(1:3)]), alpha)
      surface.list[[ln]][[i]]<- as.mesh3d(temp.ashape )
      #print(class(surface.list[[ln]][[i]]))
    }
    
  }
  return(surface.list)
}

def.surface.whole.LN <- function(cell_coords, alpha = 200){
  surface.list <- list()
  extracted.list<- list()
  idx.list <- names(cell_coords)
  for(ln in idx.list){
    print(ln)
    surface.list[[ln]] 
    extracted.list[[ln]] <- data.frame(cell_coords[[ln]])
    temp.ashape <- ashape3d(as.matrix(extracted.list[[ln]][,c(1:3)]), alpha)
    surface.list[[ln]] <- as.mesh3d(temp.ashape )
      #print(class(surface.list[[ln]][[i]]))
  }
  return(surface.list)
}



##Edge corrected volume

vol.edge.corrected <- function(surf.list, refs.list, shells){
  idx.list <- names(refs.list)
  j = 0
  temp.vol.mat <- matrix(0, nrow = sum(unlist(lapply( refs.list,nrow))), ncol = length(shells))
  for(ln in idx.list){
    print(ln)
    for(i in 1:nrow(refs.list[[ln]])){
      print(i)
      temp.vol<- apply(X = t(shells + shells[2]-shells[1]),
                       MARGIN = 2 , 
                       FUN = function(r){
                         print("a")
                         sphere.mesh <- const.sphere.mesh(r, refs.list[[ln]][i,1:3],subdivision = 4) 
                         print("b")
                         intersect.mesh <- wrap_mesh_volume(wrap_mesh_boolean(sphere.mesh,  surf.list[[ln]][[i]],op_type = "intersect"))
                         print(intersect.mesh)
                         return((intersect.mesh))
                       })
      
      j = j+1
      temp.vol.mat[j,] <- temp.vol
      
      #saveRDS(vol.surf.gast.ts , file = "surf.vol.gast.ts.Rds")
    }
  }
  return(temp.vol.mat)
}

vol.edge.corrected.whole.LN <- function(surf.list, refs.list, shells=NULL, distance = NULL){
  idx.list <- names(refs.list)
  j = 0
  
  
  if(!is.null(shells)){
    temp.vol.mat <- matrix(0, nrow = sum(unlist(lapply( refs.list,nrow))), ncol = length(shells))  
  }
  
  if(!is.null(distance)){
    temp.vol.mat <- matrix(0, nrow = sum(unlist(lapply( refs.list,nrow))), ncol = 1)  
  }
  
  for(ln in idx.list){
    print(ln)
    for(i in 1:nrow(refs.list[[ln]])){
      print(i)
      if(!is.null(shells)){
        if(1){
          #define ref cell-centered vol
          r.temp <- (shells + shells[2]-shells[1])[length(shells)]
          sphere.mesh.temp <- const.sphere.mesh(r.temp, refs.list[[ln]][i,1:3],subdivision = 4)
          surface.mesh.ref <- wrap_mesh_boolean(sphere.mesh.temp,  surf.list[[ln]],op_type = "intersect")
          
        }
        
        temp.vol<- apply(X = t(shells + shells[2]-shells[1]),
                         MARGIN = 2 , 
                         FUN = function(r){
                           #print("a")
                           sphere.mesh <- const.sphere.mesh(r, refs.list[[ln]][i,1:3],subdivision = 4) 
                           #print("b")
                           #intersect.mesh <- wrap_mesh_volume(wrap_mesh_boolean(sphere.mesh,  surf.list[[ln]],op_type = "intersect"))
                           intersect.mesh <- wrap_mesh_volume(wrap_mesh_boolean(sphere.mesh,   surface.mesh.ref,op_type = "intersect"))
                           
                           cat(r)
                           cat(" ")
                           #print(intersect.mesh)
                           return((intersect.mesh))
                         })
        cat("\n")
        
        j = j+1
        temp.vol.mat[j,] <- temp.vol
        
      }
      
      if(!is.null(distance)){
        
        #define ref cell-centered vol
        r.temp <- distance
        #print("a")
        
        sphere.mesh.temp <- const.sphere.mesh(r.temp, refs.list[[ln]][i,1:3],subdivision = 4)
        #print("b")
        temp.vol <- wrap_mesh_volume(wrap_mesh_boolean(sphere.mesh.temp,  surf.list[[ln]],op_type = "intersect"))
        #print("c")
        
        cat("\n")
        
        j = j+1
        temp.vol.mat[j,] <- temp.vol
        
      }

      
      #saveRDS(vol.surf.gast.ts , file = "surf.vol.gast.ts.Rds")
    }
  }
  return(temp.vol.mat)
}



shell.density.moving.width <- function(cell.num.per.shell, volume.per.shell, shell.range ,width){
  #shell.range <- c(1,90)
  #width <- 5
  #cell.num.per.shell <- tregs.pSTAT5.T102S.wt.num.1um.all[,-1]
  #volume.per.shell <- vol.surf.T102S.wt.delt
  vol.shell.mov.width <- matrix(0, nrow = nrow(cell.num.per.shell), ncol = length(shell.range[1]:(shell.range[2]-width+1)))
  cell.num.mov.width <- matrix(0, nrow = nrow(cell.num.per.shell), ncol = length(shell.range[1]:(shell.range[2]-width+1)))
 
  for(i in shell.range[1]:(shell.range[2]-width+1)){
    vol.shell.mov.width[,i] <- rowSums( volume.per.shell[,i:(i+width-1)])
    cell.num.mov.width[,i] <- rowSums(cell.num.per.shell[,i:(i+width-1)])
    
  }
  
  shell.density.moving.width <- cell.num.mov.width/vol.shell.mov.width
  colnames(shell.density.moving.width) <- seq(width/2,shell.range[2]-width/2,1)
  return(shell.density.moving.width)
}


avg.marker.intensity.moving.width <- function(sum.marker.per.shell, cell.num.per.shell, shell.range, width){
  #shell.range <- c(1,90)
  #width <- 5
  #cell.num.per.shell <- tregs.pSTAT5.T102S.wt.num.1um.all[,-1]
  #volume.per.shell <- vol.surf.T102S.wt.delt
  
  sum.marker.intensity.mov.width <- matrix(0, nrow = nrow(cell.num.per.shell), ncol = length(shell.range[1]:(shell.range[2]-width+1)))
  cell.num.mov.width <- matrix(0, nrow = nrow(cell.num.per.shell), ncol = length(shell.range[1]:(shell.range[2]-width+1)))
  
  for(i in shell.range[1]:(shell.range[2]-width+1)){
    sum.marker.intensity.mov.width[,i] <- rowSums( sum.marker.per.shell[,i:(i+width-1)])
    cell.num.mov.width[,i] <- rowSums(cell.num.per.shell[,i:(i+width-1)])
    
  }
  
  avg.marker.intensity.moving.width <- sum.marker.intensity.mov.width/cell.num.mov.width
  colnames( avg.marker.intensity.moving.width) <- seq(width/2,shell.range[2]-width/2,1)
  return( avg.marker.intensity.moving.width)
}

#Plotting functions
library(viridis)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


density.viridis.plot <- function(ggdata, log = F){
  
  xname <- names(ggdata)[1]
  yname <- names(ggdata)[2]
  ggdata <- ggdata[!is.na(ggdata[,2]),]
  ggdata <- cbind(ggdata, density = get_density(ggdata[,1],ggdata[,2], n = 100))
  if(!log){
    ggplot(ggdata, aes_string(xname, yname, col = "density"))+  geom_point(alpha=0.5, size = 0.3) + 
      scale_color_viridis()# +geom_abline(slope = -1.2, intercept = 1.3)
  }else{
    ggplot(ggdata, aes_string(xname, yname, col = "density"))+  geom_point(alpha=0.5, size = 0.3) + 
      scale_color_viridis()+ scale_x_continuous(trans='log10') +
      scale_y_continuous(trans='log10')
    
  }
}


marker.intensity.within.dist <- function(markers, cell.extract, cell.type, num.ref.cells ){
  center.idx <- unique(cell.extract$center.idx)
  #print(center.idx)
  temp.mat <- matrix(NA, ncol = length(markers)*2, nrow = num.ref.cells)
  colnames(temp.mat) <- rep(NA, length(markers)*2)
  counter = 0
  for(idx.markers in markers){
    counter = counter +1
    print(idx.markers)
    # 
    # data.tconv[,markers.max] <- NA
    markers.max <- paste0(cell.type, ".", idx.markers, ".max")
    markers.avg <- paste0(cell.type, ".", idx.markers, ".avg")
    # data.tconv[,markers.avg] <- NA
    # 
    #print( colnames(temp.mat)[c(counter*2-1, counter*2)])
    #print(c(markers.avg, markers.max))
    colnames(temp.mat)[c(counter*2-1, counter*2)] <- c(markers.avg, markers.max)
    #print(c(markers.avg, markers.max))
    
    temp.mat[,markers.avg][center.idx] <- sapply(center.idx,function(x){
      #print(x)
      #x = 2
      return(mean(cell.extract[,idx.markers][cell.extract$center.idx==x]))
    })
    
    temp.mat[,markers.max][center.idx] <- sapply(center.idx,function(x){
      #print(x)
      #x = 2
      return(max(cell.extract[,idx.markers][cell.extract$center.idx==x],na.rm = T))
    })
    
    # for(ref.idx in center.idx){
    #   data.tconv[,markers.max][ref.idx] <- max(tregs.extract.dist.20$whole[,idx.markers][tregs.extract.dist.20$whole$center.idx==ref.idx],na.rm = T)
    #   data.tconv[,markers.avg][ref.idx] <- mean(tregs.extract.dist.20$whole[,idx.markers][tregs.extract.dist.20$whole$center.idx==ref.idx])
    #   #print( data.tconv[,markers.avg][ref.idx])
    # }
  }
  return(temp.mat)
}

