#' Count unique trust triangles in a directed graph
#'
#' @param Adj Adjacency matrix.
#'
#' @return Integer count of trust triangles.
#' @export
Find_Number_Trust_Triangles_Unique <- function(Adj) {
  A <- sign(Adj)  # convert to 0/1
  Wedge_matrix <- A %*% A
  # Triangles are positions where both A and A^2 have edges
  triangle_matrix <- (A + Wedge_matrix) != 0 & (A == Wedge_matrix)
  NoTriangles <- sum(Wedge_matrix[triangle_matrix])
  return(as.integer(NoTriangles))
}

#' Count triangles allowing multi-edges
#'
#' @param Adj Adjacency matrix.
#'
#' @return Integer count of triangles.
#' @export
Find_Number_Triangles <- function(Adj) {
  triangle_matrix <- Adj %*% Adj %*% Adj
  No_triangles <- sum(diag(triangle_matrix)) / 3
  return(as.integer(No_triangles))
}

#' Count unique triangles in a directed graph
#'
#' @param Adj Adjacency matrix.
#'
#' @return Integer count of unique triangles.
#' @export
Find_Number_Triangles_Unique <- function(Adj) {
  A <- sign(Adj)
  triangle_matrix <- A %*% A %*% A
  No_triangles <- sum(diag(triangle_matrix)) / 3
  return(as.integer(No_triangles))
}

#' Count reciprocal 2-loops allowing multi-edges
#'
#' @param Adj Adjacency matrix.
#'
#' @return Integer count of 2-loops.
#' @export
Find_Number_2Loops <- function(Adj) {
  loop_matrix <- Adj %*% Adj
  No_loops <- sum(diag(loop_matrix)) / 2
  return(as.integer(No_loops))
}

#' Count unique reciprocal 2-loops
#'
#' @param Adj Adjacency matrix.
#'
#' @return Integer count of unique 2-loops.
#' @export
Find_Number_2Loops_Unique <- function(Adj) {
  A <- sign(Adj)
  loop_matrix <- A %*% A
  No_loops <- sum(diag(loop_matrix)) / 2
  return(as.integer(No_loops))
}

#' Count wedges allowing multi-edges
#'
#' @param Adj Adjacency matrix.
#'
#' @return Integer count of wedges.
#' @export
Find_Number_Wedges <- function(Adj) {
  wedge_matrix <- Adj %*% Adj
  No_wedges <- sum(wedge_matrix) - sum(diag(wedge_matrix))
  return(as.integer(No_wedges))
}

#' Count unique wedges in a directed graph
#'
#' @param Adj Adjacency matrix.
#'
#' @return Integer count of unique wedges.
#' @export
Find_Number_Wedges_Unique <- function(Adj) {
  A <- sign(Adj)
  wedge_matrix <- A %*% A
  No_wedges <- sum(wedge_matrix) - sum(diag(wedge_matrix))
  return(as.integer(No_wedges))
}

#' Enumerate outward trust triangles
#'
#' @param Adj Adjacency matrix.
#'
#' @return A list with the triangle count and the vertex triplets found.
#' @export
Trust_Triangles <- function(Adj) {
  # Initialize counter for number of triangles
  NoTriangles <- 0
  
  # Initialize list to store triangle triplets
  Triangle_list <- list()
  
  # Number of vertices in the graph
  n <- nrow(Adj)
  
  # Loop over each vertex v (potential "source" node)
  for (v in 1:n) {
    
    # Find all outgoing neighbors of v (v -> *)
    neigh_v <- which(Adj[v, ] != 0)
    
    # For each neighbor w of v (edge v -> w)
    for (w in neigh_v) {
      
      # Find all outgoing neighbors of w (w -> *)
      neigh_w <- which(Adj[w, ] != 0)
      
      # Find common neighbors between v and w
      # These are nodes u such that:
      # v -> u AND w -> u
      intersec <- intersect(neigh_v, neigh_w)
      
      # Each common neighbor u forms a "trust triangle":
      # v -> w, v -> u, and w -> u
      NoTriangles <- NoTriangles + length(intersec)
      
      # Store each triangle found
      if (length(intersec) > 0) {
        Triangle_list <- c(
          Triangle_list,
          lapply(intersec, function(u) c(u, v, w))
        )
      }
    }
  }
  
  # Return:
  # - total number of triangles
  # - matrix where each row is a triangle [u, v, w]
  return(list(
    NoTriangles = NoTriangles,
    Triangle_list = do.call(rbind, Triangle_list)
  ))
}

#' Enumerate directed cycle triangles
#'
#' @param Adj Adjacency matrix.
#'
#' @return A list with the cycle-triangle count and the vertex triplets found.
#' @export
Cycle_Triangles <- function(Adj) {
  # Initialize counter for number of cycle triangles
  NoTriangles <- 0
  
  # Initialize list to store triangle triplets
  Triangle_list <- list()
  
  # Transpose of adjacency matrix
  # This allows us to easily find incoming neighbors (u → v)
  AdjT <- t(Adj)
  
  # Number of vertices in the graph
  n <- nrow(Adj)
  
  # Loop over each vertex v (potential starting point of the cycle)
  for (v in 1:n) {
    
    # Outgoing neighbors of v: nodes w such that v → w
    neigh_v <- which(Adj[v, ] != 0)
    
    # Incoming neighbors of v: nodes u such that u → v
    backneigh_v <- which(AdjT[v, ] != 0)
    
    # For each outgoing neighbor w (edge v → w)
    for (w in neigh_v) {
      
      # Outgoing neighbors of w: nodes u such that w → u
      neigh_w <- which(Adj[w, ] != 0)
      
      # Find nodes u that satisfy BOTH:
      #   u → v  (from backneigh_v)
      #   w → u  (from neigh_w)
      # This completes the cycle: v → w → u → v
      intersec <- intersect(backneigh_v, neigh_w)
      
      # Each such u forms a cycle triangle
      NoTriangles <- NoTriangles + length(intersec)
      
      # Store each triangle found as (u, v, w)
      # representing the cycle: v → w → u → v
      if (length(intersec) > 0) {
        Triangle_list <- c(
          Triangle_list,
          lapply(intersec, function(u) c(u, v, w))
        )
      }
    }
  }
  
  # Each triangle is counted 3 times (once from each vertex),
  # so divide by 3 to get the correct number of unique triangles
  # 1 → 2 → 3 → 1
  # v = 1: (1,2,3)
  # v = 2: (2,3,1)
  # v = 3: (3,1,2)
  return(list(
    NoTriangles = round(NoTriangles / 3),
    Triangle_list = do.call(rbind, Triangle_list)
  ))
}

#' Enumerate wedges in a directed graph
#'
#' @param Adj Adjacency matrix.
#'
#' @return A list with the wedge count and the vertex triplets found.
#' @export
Wedges <- function(Adj) {
  # Initialize the count of wedges and a list to store them
  NoWedges <- 0
  Wedge_list <- list()
  
  n <- nrow(Adj)  # Number of vertices in the graph
  
  # Loop over each vertex v
  for (v in 1:n) {
    # Find all neighbors of v (vertices w such that v -> w exists)
    neigh_v <- which(Adj[v, ] != 0)
    
    # Loop over each neighbor w of v
    for (w in neigh_v) {
      # Find neighbors of w, excluding v to avoid trivial wedge loops
      neigh_w <- setdiff(which(Adj[w, ] != 0), v)
      
      # Loop over each neighbor u of w
      for (u in neigh_w) {
        # Each triple (v -> w -> u) forms a wedge
        NoWedges <- NoWedges + 1
        
        # Store the wedge as a vector (v, w, u)
        Wedge_list <- c(Wedge_list, list(c(v, w, u)))
      }
    }
  }
  
  # Return the total count and the full list as a matrix
  return(list(NoWedges = NoWedges, Wedge_list = do.call(rbind, Wedge_list)))
}

#' Breadth-first search over a graph
#'
#' @param Adj Adjacency matrix.
#' @param root Starting vertex index.
#'
#' @return Integer vector of visited vertices.
#' @keywords internal
BFS <- function(Adj, root) {
  N <- nrow(Adj)
  seen <- rep(FALSE, N)
  Active <- root
  Explored <- c()
  seen[root] <- TRUE
  
  while (length(Active) > 0) {
    curr <- Active[1]
    Active <- Active[-1]
    
    neigh <- which(Adj[curr, ] != 0)
    for (w in neigh) {
      if (!seen[w]) {
        Active <- c(Active, w)
        seen[w] <- TRUE
      }
    }
    
    Explored <- c(Explored, curr)
  }
  
  return(Explored)
}

#' Find strongly connected components with Tarjan's algorithm
#'
#' @param Adj Adjacency matrix.
#'
#' @return A list of strongly connected components.
#' @keywords internal
TarjanIterative <- function(Adj) {
  N <- nrow(Adj)
  index <- rep(NA, N)
  lowlink <- rep(NA, N)
  onstack <- rep(FALSE, N)
  stack <- c()
  groups <- list()
  new_index <- 1
  nextgroup <- 1
  groupid <- integer(N)
  
  sconnect <- function(v) {
    work <- list(list(v=v, i=0))
    
    while (length(work) > 0) {
      top <- work[[length(work)]]
      work <- work[-length(work)]
      v <- top$v
      i <- top$i
      
      if (i == 0) { # First visit
        index[v] <<- new_index
        lowlink[v] <<- new_index
        new_index <<- new_index + 1
        stack <<- c(stack, v)
        onstack[v] <<- TRUE
      }
      
      recurse <- FALSE
      neigh <- which(Adj[v, ] != 0)
      for (j in seq(i + 1, length(neigh))) {
        w <- neigh[j]
        if (is.na(index[w])) {
          work <- c(work, list(list(v=v, i=j)))
          work <- c(work, list(list(v=w, i=0)))
          recurse <- TRUE
          break
        } else if (onstack[w]) {
          lowlink[v] <<- min(lowlink[v], index[w])
        }
      }
      if (recurse) next
      
      if (index[v] == lowlink[v]) {
        com <- c()
        repeat {
          w <- stack[length(stack)]
          stack <- stack[-length(stack)]
          onstack[w] <<- FALSE
          com <- c(com, w)
          groupid[w] <<- nextgroup
          if (w == v) break
        }
        groups[[nextgroup]] <<- com
        nextgroup <<- nextgroup + 1
      }
      
      if (length(work) > 0) {
        w <- v
        v <- work[[length(work)]]$v
        lowlink[v] <<- min(lowlink[v], lowlink[w])
      }
    }
  }
  
  for (v in 1:N) {
    if (is.na(index[v])) sconnect(v)
  }
  
  return(groups)
}

#' Extract the giant strongly connected component
#'
#' @param Adj Adjacency matrix.
#'
#' @return Integer vector of vertex indices in the GSCC.
#' @export
GSCC <- function(Adj) {
  SCCs <- TarjanIterative(Adj)
  sizes <- sapply(SCCs, length)
  giant_idx <- which.max(sizes)
  return(SCCs[[giant_idx]])
}

#' Compute the OUT component of a directed graph
#'
#' @param Adj Adjacency matrix.
#'
#' @return Integer vector of vertices in the OUT component.
#' @export
OUT <- function(Adj) {
  giant <- GSCC(Adj)
  giant_out <- BFS(Adj, giant[1])
  Out_component <- setdiff(giant_out, giant)
  return(Out_component)
}

#' Compute the IN component of a directed graph
#'
#' @param Adj Adjacency matrix.
#'
#' @return Integer vector of vertices in the IN component.
#' @export
IN <- function(Adj) {
  trans <- t(Adj)
  giant <- GSCC(trans)
  giant_in <- BFS(trans, giant[1])
  In_component <- setdiff(giant_in, giant)
  return(In_component)
}