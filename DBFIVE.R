####################################################################################################
#' # DBFIVE: Dynamic Forwards Inference Variable Elimination for DBNs with dynamic CPTs
####################################################################################################

IDTSCALE = 0.1 # Used to generate unique id for a node = n + IDTSCALE*dt
ITSCALE = 0.001 # Used to generate unique id for a node = n + ITSCALE * t (t=0..100)
IDNSCALE = 0.01 # Used to generate unique id for a node = IDNSCALE*n + dt
MAXMEM = 16*2^30/2^20 # 16GB
####################################################################################################
#' # DBN to DBFIVE Functions
####################################################################################################
# --------------------------------------------------------------------------------------------------
#' @title getparinds
#' @description Transform vector of parent node names into parent node inds
#' @param parnamevect = char vector of par names, names include a [t-x] to indicate link to previous
#' t-slice
#' @param nnames = char vector of DBN node names
#' @return Returns parnname, parnind, parndt = vectors of par node name with [t-x] removed, parnind 
#' = index to nnames of parent node, and parndt = x
#' @export
#' @examples none
#' @name getparinds
getparinds = function(parnamevect,nnames) {
  parnname = vector('character',length(parnamevect))
  parnind = vector('integer',length(parnamevect))
  parndt = vector('integer',length(parnamevect))
  
  if(length(parnamevect)==0) {
    return(list(parnname=parnname,parnind=parnind,parndt=parndt))
  }
  for(i in 1:length(parnamevect)) {
    parname = gsub('([A-za-z0-9_ ]+)*\\[.+', '\\1',parnamevect[i]) # Extract name (inc. whitespace)
    parname = gsub('^\\s+|\\s+$','',parname) # Trim any leading/trailing whitespace
    parname = gsub(' ','_',parname) # Replace white space with _
    dtslice = gsub('.+\\[.+?([0-9]+)\\]?.+$', '\\1', parnamevect[i])#Get par node t-slice(1, 2 etc.)
    dtslice = suppressWarnings(as.integer(dtslice))
    dtslice[is.na(dtslice)] = as.integer(0) # If at same t-slice, set as 0
    ind = match(parname,nnames)
    if(!is.na(ind)) {
      parnname[i] = parname
      parnind[i] = ind
      parndt[i] = dtslice
    }
  } # end for
  return(list(parnname=parnname,parnind=parnind,parndt=parndt))
}

#' @title getchildinds
#' @description Identify child nodes for given param:
#' @param node (int) node id
#' @param maxdt (int) note that maxdt is the max number of time slices left in simulation, e.g. if 
#' only 1 slice left, any arcs to t+1 are discarded as they are outside the simulation window.
#' @param parnodeind (parnodeind (vector list of int vectors)
#' @param parnodetslice (same format as parnodeind)
#' @return Returns vector of childnodeind and childnodetslice
#' @export
#' @examples none
#' @name getchildinds
getchildinds = function(node,maxdt,parnodeind,parnodetslice) {
  xx = mapply(function(X,Y) {
    sapply(node, FUN=function(n,t) X==n & Y<=t,t=maxdt) # Up to and including maxdt t-slice
  },X=parnodeind,Y=parnodetslice)
  mask = sapply(xx,function(xx) any(xx==TRUE))
  
  ind = which(mask,arr.ind=TRUE)
  if(!is.null(ind) && length(ind)>0) {
    childnodeind = vector('integer',length(ind))
    childnodetslice = childnodeind
    for(i in 1:length(ind)) { # for each identified child node record node ind and dtslice
      childnodeind[i] = ind[i]
      childnodetslice[i] = parnodetslice[[ind[i]]][xx[[ind[i]]]]
    }
    return(list(childnodeind=childnodeind,childnodetslice=childnodetslice))
  } else {
    return(list(childnodeind=vector('integer',0),childnodetslice=vector('integer',0)))
  } #end if
}

#' @title dbn2DFnodes
#' @description Convert dbn$node name, states, parents to DBFIVE$nnames, statenames, parnodeind,
#' parnodetslice, childnodeind, childnodetslice
#' @param dbn = DBNSMILER object, 
#' @param maxdt = (int) max# dt-slices (e.g. Markov network, maxdt = 1; BN, maxdt=0)
#' @return Returns DBFIVE$nnames, statenames, parnodeind, parnodetslice, childnodeind, 
#' childnodetslice
#' @export
#' @examples none
#' @name dbn2DFnodes
dbn2DFnodes = function(dbn,maxdt) {
  nnames = vector('character',length(dbn$node))
  statenames = vector('list',length(dbn$node))
  numstates = vector('integer',length(dbn$node))
  parnodeind = vector('list',length(dbn$node)*maxdt)
  dim(parnodeind) = c(maxdt,length(dbn$node))
  parnodetslice = parnodeind
  childnodeind = parnodeind
  childnodetslice = parnodeind
  cpt = parnodeind
  
  for(i in 1:length(dbn$node))  { # Get nnames and statenames and numstates
    nnames[i] = dbn$node[[i]]$name
    statenames[[i]] = dbn$node[[i]]$states
    numstates[i] = length(statenames[[i]])
  }
  
  for(i in 1:length(dbn$node))  { # Get cpt foreach node
    for(j in 1:maxdt) { # foreach dt-slice
      k = min(length(dbn$node[[i]]$cpt),j)
      cpt[[j,i]] = dbn$node[[i]]$cpt[[k]]
    }
  }
  
  for(i in 1:length(dbn$node))  { # Get parnodeind foreach node
    for(j in 1:maxdt) { # foreach dt-slice
      k = min(length(dbn$node[[i]]$parnodenames),j)
      x=getparinds(dbn$node[[i]]$parnodenames[[k]],nnames)
      parnodeind[[j,i]] = x$parnind
      parnodetslice[[j,i]] = x$parndt
    }
  }
  
  for(i in 1:length(dbn$node)) { # Get childnodeind foreach node
    for(j in 1:maxdt) { # foreach dt-slice
      x=getchildinds(i,j-1,parnodeind[j,],parnodetslice[j,])
      childnodeind[[j,i]] = x$childnodeind
      childnodetslice[[j,i]] = x$childnodetslice
    } # end for times t
  }# end foreach node
  
  return(list(nnames=nnames,statenames=statenames,parnodeind=parnodeind,parnodetslice=parnodetslice,
            childnodeind=childnodeind,childnodetslice=childnodetslice,numstates=numstates,cpt=cpt))
}

#...................................................................................................
# # Test
# x = dbn2DFnodes(dbn,2)
# nnames = x$nnames
# statenames = x$statenames
# parnodeind = x$parnodeind
# parnodetslice = x$parnodetslice
# childnodeind = x$childnodeind
# childnodetslice = x$childnodetslice
# numstates = x$numstates
# cpt = x$cpt


# --------------------------------------------------------------------------------------------------
#' @title dbn2DFevid
#' @description Convert dbn$evid structure to DBFIVE$evid_t structure and also adds 'Time_of_Year'
#' evidence
#' @param dbn DBNSMILER object
#' @param dbn$evid DBNSMILER evid object
#' @param nnames DBFIVE char vec of node names
#' @return DBFIVE$evid_t[[t,node]][states] structure
#' @export
#' @examples none
#' @name dbn2DFevid
dbn2DFevid = function(dbn,E,nnames) {
  globdater = dbn$globdater
  nt = dbn$globdater$num[2] - dbn$globdater$num[1] + 1 # Total #t-slices in simulation window
  evid_t = vector('list',length(dbn$node)*nt)
  dim(evid_t) = c(nt,length(dbn$node))
  
  # Generate and add datetime 'Time_of_Year' evidence
  chartime = num2chardate(dbn$globdater$num[1]:dbn$globdater$num[2],dbn$datetype)
  dtnode = 'Time_of_Year'
  ind = findnodebyname(dbn,dtnode)
  if(ind>0) {
    E=gendatetimeevid(E,dtnode,chartime,1,dbn$node[[ind]]$states,dbn$datetype)
  } else { # Otherwise do nothing as there is no dtnode
  }
  if(nrow(E)<=0)
    return(evid_t)
  for(i in 1:nrow(E)) {
    nodeind = match(E$Node[i],nnames)
    if(length(nodeind)==1 && !is.na(nodeind)) { # Check evidence node is in the DBN, stop if not
      # First convert evidence into a probability beliefs vector
      beliefvect = vector('numeric',length(dbn$node[[nodeind]]$states))
      if(E$SelectState[i]!='NA') { # Set state evidence
        beliefvect[match(E$SelectState[i],dbn$node[[nodeind]]$states)] = 1
      } else {
        beliefvect = as.numeric(unlist(strsplit(E$SoftEvidence[i],split=',')))
      }
      
      # Handle time slicing
      t = chardate2num(c(E$DateStart[i],E$DateEnd[i]), dbn$datetype)
      
      # Intersect trange with globdater trange - 3 cases
      if(t[1]<globdater$num[1] && t[2]>=globdater$num[1] && t[2]<=globdater$num[2]) {
        trange = c(globdater$num[1],t[2])
      } else if(t[1]>=globdater$num[1] && t[1]<=globdater$num[2] &&
                  t[2]>=globdater$num[1] && t[2]<=globdater$num[2]) {
        trange = c(t[1],t[2])
      } else if(t[1]>=globdater$num[1] && t[1]<=globdater$num[2] && t[2]>globdater$num[2]) {
        trange = c(t[1],globdater$num[2])
      } else {
        trange = NA
      }
      
      if(!any(is.na(trange))) {
        # Convert trange to t-slice range (which goes from 1 to numslices)
        trangevect = (trange[1]-globdater$num[1]+1):(trange[2]-globdater$num[1]+1)
        #       cat('Evid node',E$Node[i],trangevect,'\n')
        for(j in 1:length(trangevect)) { # Set the evidence!
          # Note that trangevect[j]= t-slice, nodeind = indmatch col as nnames in same order as 
          # dbn$node
          evid_t[[trangevect[j],nodeind]] = beliefvect
          #cat('Node:',nnamesur[ind],':',beliefcpt,'\n') # DEBUGGING
        } # end for each t-slice in trangevect
      }
    } else {
      cat('Node:',E$Node[i],', Evidence node not found in DBN\n')
    } # end if evid valid
  } # end foreach row of evidence
  return(evid_t)
}

#...................................................................................................
# # Test
# evid_t = dbn2DFevid(dbn,E,nnames)



####################################################################################################
#' # Inferencing functions
####################################################################################################

# --------------------------------------------------------------------------------------------------
#' @title calcnodecost
#' @description Calculates cost of doing variable elimination for a given node - based on max dist 
#' from root nodes. To help differentiate between same distance nodes, add to the cost the 
#' computational cost. So cost = <distance as integer, 0=root node>.<computational cost on 
#' [0.1,0.9]>. This means all parents of node A have a cost such that floor(a)>floor(parent_node(a))
#' <br/>. Let b<-a denote b as a child of a, if b<-a and c<-a but b<-c, then dist(a,b,c) = (0,2,1)
#' <br/>ASSUMES: DAGs, so max dist is finite
#' @param parnodeind ([[1..maxdt,node]][parent node ind to nnames])
#' @param parnodetslice ([[1..maxdt,node]][par k for $\Pi(t-k)$])
#' @param childnodeind ([[1..maxdt,node]][child node ind to nnames])
#' @param childnodetslice ([[1..maxdt,node]][child k for $\Pi(t+k)$])
#' @param numstates (vec of number of node states)
#' @param maxdt = (int) max# dt-slices (e.g. Markov network, maxdt = 1; BN, maxdt=0) 
#' @return list of $nodedist[node] and $nodecost[node] as vectors
#' @export
#' @examples none
#' @name calcnodecost
calcnodecost = function(parnodeind,parnodetslice,childnodeind,childnodetslice,numstates,maxdt) {
  # 1) Calculate nodecost based on estimated size of CPD at a node: cost is on [0.1,0.9]
  nc = matrix(0,nrow=maxdt,ncol=length(numstates))
  for(k in 1:maxdt) {
    for(i in 1:ncol(parnodeind)) {
      # Calculated expanded joint probability CPD size:
      parcost = 1
      if(length(parnodeind[[maxdt,i]])>0) {
        parcost = prod(numstates[parnodeind[[maxdt,i]]]) 
      }
      childcost = 1
      if(length(childnodeind[[k,i]])>0) {
        childcost = prod(
          sapply(childnodeind[[k,i]],function(X,P) prod(numstates[P[[X]]]), P=parnodeind[maxdt,]))
        childcost = childcost / (numstates[i]^length(childnodeind[[k,i]]))
      }
      nc[k,i] = numstates[i] * parcost * childcost
    }
  }
  nc = 0.1 + log(nc)/max(log(nc)) * 0.8
  
  # 2) Calculate nodecost for nodes at a given time slice
  # Note: for nodes which have at least 1 child at a different t-slice t+m, all effects of these 
  # arcs can be ignored iff we eliminate all nodes in all t-slices t+m before eliminating nodes in
  # t. Therefore, we treat all nodes with no children (leaf nodes) and all nodes with children at 
  # t-slice t+m (terminal leaf nodes at t) as having a node distance nd of 0.
  nd = matrix(Inf,nrow=maxdt,ncol=length(numstates))
  for(k in 1:maxdt) {
    ndprev = -nd[k,] # previous dist matrix
    d = 0
    ind = which(sapply(childnodeind[k,],FUN=function(X) length(X)==0),arr.ind=TRUE)
    ind = c(ind, which(sapply(childnodetslice[k,],FUN=function(X)any(X>0)),arr.ind=TRUE))
    while(!all(nd[k,]==ndprev)) {
      ndprev = nd[k,]
      nd[k,ind] = d # Set distance from terminal leaf node
      #Get vector of parnodes ofall current nodes in ind, filter to only keep those in same t-slice
      ind = unlist(sapply(ind,function(X,Y,Z) Y[[X]][Z[[X]]==0],
                          Y=parnodeind[k,],Z=parnodetslice[k,]))
      d = d + 1
    }
  }
  nodecost = nc + nd

  # 3) Calculate distance from root nodes - equiv to (2) above but in other direction
  nd2 = matrix(Inf,nrow=maxdt,ncol=length(numstates))
  for(k in 1:maxdt) {
    ndprev = -nd2[k,] # previous dist matrix
    d = 0
    ind = which(sapply(parnodeind[k,],FUN=function(X) length(X)==0),arr.ind=TRUE)
    ind = c(ind, which(sapply(parnodetslice[k,],FUN=function(X)any(X>0)),arr.ind=TRUE))
    while(!all(nd2[k,]==ndprev)) {
      ndprev = nd2[k,]
      nd2[k,ind] = d # Set distance from root node
      #Get vector of parnodes ofall current nodes in ind, filter to only keep those in same t-slice
      ind = unlist(sapply(ind,function(X,Y,Z) Y[[X]][Z[[X]]==0],
                          Y=childnodeind[k,],Z=childnodetslice[k,]))
      d = d + 1
    }
  }

  return(list(nodecost=nodecost,nodedist=nd,nodedist2=nd2))
}
#...................................................................................................
# # Test
# x = calcnodecost(parnodeind,parnodetslice,childnodeind,childnodetslice,numstates,2)
# nodecost = x$nodecost
# nodeorder = x$nodeorder

# --------------------------------------------------------------------------------------------------
#' @title getparpars
#' @description Create parpars structure of node,parents, grandparents, great- etc. foreach node
#' @param parnodeind ([[1..maxdt,node]][parent node ind to nnames])
#' @param parnodetslice ([[1..maxdt,node]][par k for $\Pi(t-k)$])
#' @param maxdt = (int) max# dt-slices (e.g. Markov network, maxdt = 1; BN, maxdt=0)
#' @param numt is number of t-slices to generate parpars
#' @return parpars$n[[maxdt,numnodes]]=vector of nodeinds, $t[[maxdt,numnodes]]=vec of twindow index
#' @export
#' @examples none
#' @name getparpars
getparpars = function(parnodeind,parnodetslice,maxdt,numt) { 
  x = vector('list',numt*ncol(parnodeind))
  xi = sapply(x,function(X) X=integer(0))
  dim(xi) = c(numt,ncol(parnodeind))
  parpars = list(n=xi, t=xi)
  pp = NULL
  for(i in 1:numt) {
    for(j in 1:ncol(parnodeind)) {
      pp$n = j # Initialise for getting parents
      pp$t = i
      toffset = i # Need offset as parents are specified as 0=same t-slice, 1 = prev t-slice etc.
      # So to get t-index to time window, need to offset by previous pp$t
      while(length(pp$n)>0) {
        parpars$n[[i,j]] = c(parpars$n[[i,j]], pp$n)
        parpars$t[[i,j]] = c(parpars$t[[i,j]], pp$t)
        ppn = as.vector(unlist(sapply(1:length(pp$n),function(X,pp,P,maxdt) 
          P[[min(pp$t[X],maxdt),pp$n[X]]], pp=pp,P=parnodeind,maxdt=maxdt)))
        ppt=as.vector(unlist(sapply(1:length(pp$t),function(X,pp,P,toff,maxdt)
          toff[X]-P[[min(pp$t[X],maxdt),pp$n[X]]],pp=pp,P=parnodetslice,toff=toffset,maxdt=maxdt)))
        toffset = ppt
        pp$n = ppn
        pp$t = ppt
      }
      ppx = parpars$n[[i,j]] + IDTSCALE*parpars$t[[i,j]]
      ind = match(unique(ppx),ppx) # Remove any duplicate nodes/linked nodes
      parpars$n[[i,j]] = parpars$n[[i,j]][ind]
      parpars$t[[i,j]] = parpars$t[[i,j]][ind]
    } # end foreach node
  } # end foreach parnodetslice scenario
  return(parpars)
} 

# --------------------------------------------------------------------------------------------------
#' @title getnodeorder
#' @description Create nodeorder structure foreach node,k scenario - i.e. sequence in which to
#' eliminate nodes to get marginal for a given target node and dt-slice.
#' @param parpars$n[[k=1..maxdt for $\Pi(t-k)$,node]] vec of node inds and $t[[k,node]] of node 
#' dt-slice t-k
#' @param childnodeind ([[1..maxdt,node]][child node ind to nnames])
#' @param childnodetslice ([[1..maxdt,node]][child k for $\Pi(t+k)$])
#' @param nodecost[1..maxdt,node] containing cost value $OC=D+C$
#' @param maxdt = (int) max# dt-slices (e.g. Markov network, maxdt = 1; BN, maxdt=0)
#' @param nnames (char vec of node names)
#' @param customPi function link to custom $\Pi'$ function
#' @return nodeorder$n[[maxdt,numnodes]]=vector of nodeinds, $k[[maxdt,numnodes]]=vec of k
#' @export
#' @examples none
#' @name getnodeorder
getnodeorder = function(parpars,childnodeind,childnodetslice,nodecost,maxdt,nnames,customPi) {
  # nodeorder = $\Pi(t-k)$, nodeorder[[k\in K,target nodes]] = vec of nodes sorted in order to elim
  x = vector('list',maxdt*length(nnames))
  xi = sapply(x,function(X) X=integer(0))
  dim(xi) = c(maxdt,length(nnames))
  nodeorder = list(n=xi,k=xi)
  
  # Get parpars$n/t[[maxdt,node]], let tau = parpars$t[[maxdt,node]]; then tau=maxdt corresponds to 
  # k=0, tau=maxdt-1->k=1, tau=maxdt-2->k=2
  for(j in 1:length(nnames)) { 
    ppn = parpars$n[[maxdt,j]]
    ppt = parpars$t[[maxdt,j]]
    
    if(is.null(customPi)) {
      # Create $\Pi(t-k)$ - this works for general $\Pi(t-k)$
      for(i in 1:maxdt) { # foreach k value, k=0...K, i=k+1
        if(i==1) {
          nodeorder$n[[i,j]] = ppn[ppt==maxdt-i+1][-1] # drop first element in t-k, k=0 slice as 
        } else {                                       # that is the TARGET node! don't eliminate!
          nodeorder$n[[i,j]] = ppn[ppt==maxdt-i+1]
        }
        #nodeorder$n[[i,j]] = ppn[ppt==maxdt-i+1]
        if(length(nodeorder$n[[i,j]])>0) {
          # Re-order nodeorder by nodecost
          nc = nodecost[nodeorder$n[[i,j]]] # deliberately only use nodecost[1,j]
          nodeorder$n[[i,j]] = nodeorder$n[[i,j]][order(nc)]
          nodeorder$k[[i,j]] = vector('integer',length(nodeorder$n[[i,j]])) #vec of zeros
        }
      } # end foreach k
    } else {
      nodeorder = customPi(nodeorder,j,ppn,ppt,maxdt,childnodeind,childnodetslice,nnames)
    }
  } # end foreach node
  return(nodeorder)
}

# --------------------------------------------------------------------------------------------------
#' @title filterchildnode
#' @description Filter child nodes to create a customised list foreach node in nodeorder. When doing
#' backwards inference, child nodes can only be in $\Pi(t)$ since all child nodes at t+m have 
#' already been eliminated. For given $k$, filter ONLY by nodeind $\Pi(t-k)$
#' @param childnodeind ([[1..maxdt,node]][child node ind to nnames])
#' @param childnodetslice ([[1..maxdt,node]][child k for $\Pi(t+k)$])
#' @param nodeorder$n[[k=1..maxdt for $\Pi(t-k)$,target node]] vec of node inds in order to elim
#' and $t[[k,target node]] vec of node dt-slice t-k
#' @param maxdt = (int) max# dt-slices (e.g. Markov network, maxdt = 1; BN, maxdt=0) 
#' @return nodeordch[[maxdt,numnodes]][[foreach nodeorder node]]=vector of child nodeinds
#' @return nodeordcht ditto nodeordch but vector of child node dt-slices
#' @return cptaddch ditto but vector of nodeorder[[i,j]][l] node + child nodeinds of child CPT's to 
#' add for multiply step. Already takes into account previously eliminated t-slices AND eliminated
#' nodes in this $\Pi'$.
#' @return cptaddcht ditto cptaddch but childnode dtslice
#' Note1: that separate nodeordch and cptaddch are needed as nodeordch is needed to populate linkind
#' even if no new child cpt's are being added <br/>
#' @export
#' @examples none
#' @name filterchildnode
filterchildnode = function(childnodeind,childnodetslice,nodeorder,targetnodes,maxdt) {
  nodeordch =vector('list',maxdt*ncol(childnodeind))#nodeorder child node n foreachnode in nodeorder
  dim(nodeordch) = c(maxdt,ncol(childnodeind))
  nodeordcht = nodeordch # child node dt-slice
  cptaddch = nodeordch # cptadd = new child node cpt's to add (for the multiply step of VE)
  cptaddcht = nodeordch # same as cptaddch but with child node dtslice
  
  for(j in 1:ncol(nodeorder$n)) { #foreach node
    for(i in 1:nrow(nodeorder$n)) { #foreach k
      ppn=nodeorder$n[[i,j]] # Use ppn and ppk to get child nodes
      ppk=nodeorder$k[[i,j]] 
      tn = targetnodes[[i,j]] # Target nodes 
      # Get list of childnodes - use childnodeind/tslice[maxdt,j] to get *max* set of child nodes 
      # and subset later for child nodes that exceed tend since nodeordch[k,j] and cptaddch[k,j] are
      # used forall t=1...tend.
      C=lapply(ppn,function(X,CNI,dt) CNI[[dt,X]],CNI=childnodeind,dt=i)#childnodeinds list
      CT=lapply(ppn, function(X,dt,ct) childnodetslice[[dt,X]],dt=i,ct=childnodetslice)#chnodetslces
      if(length(C)>0)
        CNT=lapply(1:length(C),function(X,C,CT,ppk)C[[X]]+IDTSCALE*(CT[[X]]+maxdt-ppk[X]),
                   C=C,CT=CT,ppk=ppk)
      nodeordch[[i,j]] = vector('list',length(C))
      nodeordcht[[i,j]] = vector('list',length(C))
      cptaddch[[i,j]] = vector('list',length(C))
      cptaddcht[[i,j]] = vector('list',length(C))
      nodesdone = vector('numeric',0) # empty vector of node cpts already 'done' (multiplied)
      for(k in 1:i) { 
        # Coz backwards inf, pre-populate nodesdone with ppn+IDTSCALE*(maxdt+k-ppk) - already elim'd
        nodesdone = c(nodesdone,tn + IDTSCALE*(maxdt+k)) # filter out this node j at future t-slices
                                          # since already incorporated/elim'ed when infer backwards
        nodesdone = c(nodesdone,unlist(lapply(nodeorder$n[,j],function(X)X)) + 
                        IDTSCALE*(maxdt+k-unlist(lapply(nodeorder$k[,j],function(X)X)) )  )
      }
      if(i==1) { # if target t-slice t-k,k=0 - add the target node
        ppnf = c(ppn,tn) # use ppnf to filter child nodes
        ppf = c(ppn,tn)+IDTSCALE*(maxdt-c(ppk,0*tn)) # ppf is ppnf + IDTSCALE*dt
      } else { # otherwise, do not add target node
        ppnf = ppn
        ppf = ppn + IDTSCALE*(maxdt-ppk)
      }
      if(length(C)!=0) {
        for(l in 1:length(C)) {
          nodeordch[[i,j]][[l]] =  C[[l]][C[[l]] %in% ppnf]
          nodeordcht[[i,j]][[l]] = CT[[l]][C[[l]] %in% ppnf]
          mask = c(ppn[l],C[[l]]) %in% ppnf & !(c(ppf[l],CNT[[l]]) %in% nodesdone)
          cptaddch[[i,j]][[l]] = c(ppn[l],C[[l]])[mask]
          cptaddcht[[i,j]][[l]] = c(-ppk[l],CT[[l]])[mask]
          nodesdone = c(nodesdone,CNT[[l]],ppf[l]) #add to nodesdone those that have been multiplied
        } # end foreach nodeorder node
      } else {
        nodeordch[[i,j]] = vector('list',0) # initialise to empty lists
        nodeordcht[[i,j]] = vector('list',0)
        if(i==1) {
          cptaddch[[i,j]] = tn # If no child length(C)==0, add the target node anyway if i==1
          cptaddcht[[i,j]] = 0*tn
        } else {
          cptaddch[[i,j]] = vector('list',0)
          cptaddcht[[i,j]] = vector('list',0)
        }
      }
    } # end foreach k
  } # end foreach node
  return(list(nodeordch=nodeordch,nodeordcht=nodeordcht,cptaddch=cptaddch,cptaddcht=cptaddcht))
}


# --------------------------------------------------------------------------------------------------
#' @title filtercpt
#' @description load from cpt into filtcpt, filtering by evidence and propagating
#' filtered evidence (P(X=x)=0) onto child nodes.
#' @param  filtcpt[[1..maxdt,node]][numeric vector of CPT] - CPT to be filtered
#' @param nnames (char vec of node names)
#' @param maxdt (int) max# dt-slices (e.g. Markov network, maxdt = 1; BN, maxdt=0)
#' @param parnodeind ([[1..maxdt,node]][parent node ind to nnames])
#' @param parnodetslice ([[1..maxdt,node]][par k for $\Pi(t-k)$])
#' @param numstates (vec of number of node states)
#' @param childnodeind ([[1..maxdt,node]][child node ind to nnames])
#' @param childnodetslice ([[1..maxdt,node]][child k for $\Pi(t+k)$])
#' @param evid_t[[1..tend,node]][vector of len numstates containing node marginals]
#' @param cpt ([[1..maxdt,node]][numeric vector of CPT])
#' @param tend = index to last t-slice in simulation window 1:tend
#' @return Jcpt[[1..tend,target node]] with $p = numeric vector of p, $n = int vector of nodes n in
#' $J$ in order in which $J$ is composed, $t = int vector of node t-slices 1..tend, and
#' $s = list of state-labels (as indices) foreach node in $n$
#' @export
#' @examples none
#' @name filtercpt
filtercpt = function(filtcpt,nnames,maxdt,parnodeind,parnodetslice,numstates,childnodeind,
                     childnodetslice,evid_t,cpt,tend) {
  # Pre-load CPT with evidence into 'filtered cpt' filtcpt
  # 1) Populate filtcpt ignoring evidence
  for(i in 1:tend) { #foreach t-slice 1..tend
    for(j in 1:length(nnames)) { # foreach node
      kpar = min(i,maxdt) # indexes K t-slice for parnodeind and cpt
      filtcpt[[i,j]]$p = cpt[[kpar,j]]
      filtcpt[[i,j]]$n = c(parnodeind[[kpar,j]],j)
      filtcpt[[i,j]]$t = i - c(parnodetslice[[kpar,j]],0)
      filtcpt[[i,j]]$s = lapply(numstates[filtcpt[[i,j]]$n],function(X)1:X)
    }
  }
  # 2) Filter by evidence
  for(i in 1:tend) { #foreach t-slice 1..tend
    for(j in 1:length(nnames)) { # foreach node
      if(!is.null(evid_t[[i,j]])) { # if evid exists
        kchild = min(tend-i+1,maxdt) # indexes K t-slice for childnodeind
        n = c(j,childnodeind[[kchild,j]]) # vector of nodes containing evidence node
        t = i+c(0,childnodetslice[[kchild,j]]) # vector of t-slices corresponding to n
        E = evid_t[[i,j]]
        maskind = which(E==0) # index to zero elements in evidence
        
        for(l in 1:length(n)) { # foreach node in n
          # Get index njind to n[l]'s CPD for position of evid/cpt node j
          njind = which(filtcpt[[t[l],n[l]]]$n==j & filtcpt[[t[l],n[l]]]$t==t[1])#t[1]=evid node
          if(t[l]>0 && t[l]<=tend & length(njind)==1) { # ensure t is within time bounds
            ns = sapply(filtcpt[[t[l],n[l]]]$s,function(X)length(X)) # numstates vector for node n[l]
            lentrues =prod( c(ns,1)[(njind+1):max(length(ns),njind+1)] )#width of a block of trues
            jns = ns[njind] # number of states of node j in node n[l]'s CPD
            rmmask = rep(FALSE,prod(ns))
            for(m in 1:jns) { # foreach state in evidence node
              # Note, even though m may not correspond to the actual labels of the states, it is 
              # necessary to use 1:ns as labels to ensure correct length/positioning of lentrues mask
              mask = rep(c(rep(FALSE,(m-1)*lentrues),rep(TRUE,lentrues),rep(FALSE,(jns-m)*lentrues)),
                         length.out=prod(ns)) # Mask out state m in CPD
              if(l==1) {# Only for the evidence node n[l] (not child n), write evidence to this block
                if(E[m]<0)
                  filtcpt[[t[l],n[l]]]$p[mask] = -E[m]#Hard set CPT by evid if evid was -ve
                else
                  filtcpt[[t[l],n[l]]]$p[mask] = filtcpt[[t[l],n[l]]]$p[mask] * E[m]#Scale CPT by evid
                # filtcpt[[t[l],n[l]]]$p[mask] = filtcpt[[t[l],n[l]]]$p[mask] * E[m]#Scale CPT by evid
              }
              if(m %in% maskind)
                rmmask = rmmask | mask # add to removal mask
            }
            filtcpt[[t[l],n[l]]]$p = filtcpt[[t[l],n[l]]]$p[!rmmask] # Filter out evid P(X=x)=0 states
            filtcpt[[t[l],n[l]]]$s[[njind]] = setdiff(1:jns,maskind) # Write filtered numstates
          }
          # Note that $t and $n remain unchanged as actual par-child relationships and time unaffected
        } # end foreach node in n to filter
      } # end evid exists
    } # end foreach node
  } # end foreach t
  return(filtcpt)
}

# --------------------------------------------------------------------------------------------------
#' @title multiply 
#' @description multiply two Joint CPDs ($p=numeric vec, $n and $t int vec, 
#' $s=list vec of statelabels foreach node in n); equiv to outer product filtering by common
#' nodes nc to both J's.
#' Consider common nodes $nc\in n1,n2$ and CPD of n1: p1. For each state combination of nc=
#' <nc[1]=state_i1,nc[2]=state_i2,...>, a mask matching those elements in p1 contains all possible
#' state combinations of $n1 \ nc$. This occurs because of the structure of p1; starting with n[1],
#' forall elements in p1 such that n[1]=state_1, we get every combination of 
#' n[2]=state_1..numstates2, n[3]=state_1..numstate3 etc. This is repeated for n[1]=state_2,3... and
#' foreach successive node n[2], n[3]...<br/>
#' @param p1 and p2 of the two Joint CPDs to multiply of form numeric vector of p
#' @param n1 and n2 of the two Joint CPDs to multiply of form int node id in order in which J is
#' composed
#' @param t1 and t2 of the two Joint CPDs to multiply of form int t-slices 1..tend
#' @param s1 and s2 of the two Joint CPDs to multiply of form list of state-label ids
#' @return multiplied result in p,n,t,s list structure (same as input parameters and J structure)
#' @export
#' @examples none
#' @name multiply
multiply = function(p1,p2,n1,n2,t1,t2,s1,s2) {
  # 1. Restructure p,n,t,s by putting common nodes in n1, n2 into the 'bottommost level' i.e. the
  # last elements in n1/n2
  nt1 = n1 + ITSCALE*t1 # Get unique (n,t) nodes nt
  nt2 = n2 + ITSCALE*t2
  nc = intersect(nt1,nt2) # Common nodes in both n1 and n2
  ns1 = sapply(s1,function(X)length(X)) # numstates of each node in p1
  ns2 = sapply(s2,function(X)length(X)) # ditto for p2
  
  mask1 = !(nt1 %in% nc) # Prepare outputs
  mask2 = !(nt2 %in% nc)
  n = c(n1[!mask1],n1[mask1],n2[mask2]) # Outputs n, t, and s; for both nc empty or non-empty
  t = c(t1[!mask1],t1[mask1],t2[mask2])  # order is nc nodes first, then uniq nodes in n1 then n2
  s = c(s1[!mask1],s1[mask1],s2[mask2])  # this means when doing %o%, must do p2 %o% p1
  p = vector('numeric',0)
  if(length(nc)==0) { # No intersect between p1 and p2 - can use outer product directly
    p = as.vector(p2 %o% p1)
  } else { # There are common nodes nc in p1 and p2
    ind1 = sapply(rev(nc),function(X,nt1)which(X==nt1),nt1=nt1) # reverse nc coz expand.grid makes
    ind2 = sapply(rev(nc),function(X,nt2)which(X==nt2),nt2=nt2) # cols in reverse order
    lentru1=vapply(ind1,function(X,ns)prod(c(ns,1)[(X+1):max(length(ns),X+1)] ),ns=ns1,1)
    lentru2=vapply(ind2,function(X,ns)prod(c(ns,1)[(X+1):max(length(ns),X+1)] ),ns=ns2,1)
    # Create all combo of stateinds to use as offset. stateind1,2 will have same nrows 
    # as both come from nc; also re-generate state labels as 1:ns to index correctly in mask1/2
    stateind1 = as.matrix( expand.grid(lapply(ind1,function(X,ns1)1:ns1[X],ns1=ns1)) )
    stateind2 = as.matrix( expand.grid(lapply(ind2,function(X,ns2)1:ns2[X],ns2=ns2)) )

    for(i in 1:nrow(stateind1)) {                  # as both come from nc
      mask1 = rep(TRUE,length(p1))
      mask2 = rep(TRUE,length(p2))
      for(X in 1:length(nc)) { # Build up mask foreach combo of all nc nodestates
        mask1 = mask1 & rep( c(rep(FALSE,(stateind1[i,X]-1)*lentru1[X]),rep(TRUE,lentru1[X]),
          rep(FALSE,(ns1[ind1[X]]-stateind1[i,X])*lentru1[X])), length.out=length(p1))
        mask2 = mask2 & rep( c(rep(FALSE,(stateind2[i,X]-1)*lentru2[X]),rep(TRUE,lentru2[X]),
          rep(FALSE,(ns2[ind2[X]]-stateind2[i,X])*lentru2[X])), length.out=length(p2))
      }
      p = c(p,as.vector(p2[mask2] %o% p1[mask1]))
    }
  }
  return(list(p=p,n=n,t=t,s=s))
}


# --------------------------------------------------------------------------------------------------
#' @title marginalise
#' @description marginalise - marginalise over specified nodes nvect,tvect
#' @param jp numeric vector of p (from J structure)
#' @param n int vector of nodes n in $J$ in order in which $J$ is composed
#' @param t int vector of node t-slices 1..tend
#' @param s list of state-labels (as indices) foreach node in $n$
#' @param nvect vector of node ids to marginalise over
#' @param tvect vector of t-slice 1..tend to marginalise over
#' @return marginalised result in p,n,t,s list structure (same as input parameters and J structure)
#' @export
#' @examples none
#' @name marginalise
marginalise = function(jp,n,t,s,nvect,tvect) {
  for(i in 1:length(nvect)) {
    ns = vapply(s,function(X)length(X),1) #num states foreach node in Jn
    # For node nvect[i] to marginalise, build masks for positions in p-vect foreach state state_j
    # and then sum p(mask_j). Mask is based on lentrues = (#cells that are TRUE in one period) and
    ind = which(nvect[i]==n & tvect[i]==t,arr.ind=TRUE)
    lentrues = prod( c(ns,1)[(ind+1):max(length(ns),ind+1)] )
    p = vector('numeric',prod(ns)/ns[ind])
    for(j in 1:ns[ind]) {
      mask = rep( c(rep(FALSE,(j-1)*lentrues),rep(TRUE,lentrues),
                    rep(FALSE,(ns[ind]-j)*lentrues)), length.out=length(jp))
      p = p + jp[mask]
    }
    n = n[-ind]
    t = t[-ind]
    s = s[-ind]
    jp = p
  }
  return(list(p=p,n=n,t=t,s=s))
}

# --------------------------------------------------------------------------------------------------
#' @title oneslicesumprod
#' @description For a specified slice e.g. $\Pi'(t)$, find the partial marginal probability 
#' e.g. $P(\pi'(t)|\pmb{E})$ as defined by nodes to eliminate lknodeorder (and corresponding CPTs to
#' add lkcptaddch/t). 
#' @param t = current time slice
#' @param Jb joint struct to initialise backwards sweep with (NULL if nothing to initialise), 
#' $p numeric vec of p, $n int vec of nodes n in $J$ in
#' order in which $J$ is composed, $t int vec of node t-slices 1..tend, $s list of state-labels (as 
#' indices) foreach node in $n$
#' @param indb - vector of indices to Jb$n/t for nodes to marginalise
#' @param Jf - like Jb but to "initialise" (elim at end) for forwards sweep due to 
#' backwards elimination order
#' @param indf - vector of indices to Jf$n/t for nodes to marginalise
#' @param Jcpt[[1..tend,target node]] is the original joint populated by $cpt (same struct as J) 
#' @param lknodeorder = $\Pi'(t)$ defines the order of elimination of nodes in the slice in order
#' to get slice marginals e.g. $P(\pi'(t)|\pmb{E})$. lknodeorder $n = node id, $k = dt-slice
#' @param lkcptaddch = list of vec of nodeinds corresponding to CPT's to add for multiply foreach
#' elimination step in lknodeorder
#' @param lkcptaddcht = same as lkcptaddch but for node dt-slice
#' @return J joint structure with $\Pi(\pi'(t)|\pmb{E})$; $p,$n,$t,$s (same struct as Jinit)
#' @export
#' @examples none
#' @name oneslicesumprod
oneslicesumprod = function(t,Jb,indb,Jf,indf,Jcpt,lknodeorder,lkcptaddch,lkcptaddcht) {
  J = NULL
  if(!is.null(Jb)) { # Nodes to elim in Jb before getting started; e.g. $P(\pi(t+1)|E)$
    J = marginalise(Jb$p,Jb$n,Jb$t,Jb$s,Jb$n[indb],Jb$t[indb])
  }
  for(i in 1:length(lkcptaddch)) { # Eliminate nodes in order lknodeorder (equiv. lkcptaddch)
    if(length(lkcptaddch[[i]])>0) { # There are nodes to multiply
      for(j in 1:length(lkcptaddch[[i]])) { # foreach node to multiply
        if(i==1 && j==1 && is.null(J)) { # For first node in elimination
          J = Jcpt[[t+lkcptaddcht[[i]][j],lkcptaddch[[i]][j]]]
        } else { 
          tj = t+lkcptaddcht[[i]][j]
          nj = lkcptaddch[[i]][j]
          J = multiply(J$p,Jcpt[[tj,nj]]$p,J$n,Jcpt[[tj,nj]]$n,
                       J$t,Jcpt[[tj,nj]]$t,J$s,Jcpt[[tj,nj]]$s)
        }
      } # end foreach node to multiply
    } # end if there are nodes to multiply
    # cat('J size at i=',i,' is ',length(J$p),'\n')
    # Assumes there is always a node to eliminate
    if(length(lknodeorder$n)>0)
      J = marginalise(J$p,J$n,J$t,J$s,lknodeorder$n[i],t-lknodeorder$k[i]) # marginalise
  } # end foreach nodeorder node
  if(!is.null(Jf)) { # Eliminate Jf
    J = multiply(J$p,Jf$p,J$n,Jf$n,J$t,Jf$t,J$s,Jf$s)
    J = marginalise(J$p,J$n,J$t,J$s,Jf$n[indf],Jf$t[indf])
  }
  if(sum(J$p)<1)
    J$p = J$p / sum(J$p) # Normalise = P(\theta\cap Y)/P(Y)
  return(J)
}

#' @title oneslicemargs
#' @description For a specified slice e.g. $\Pi'(t)$, find the posterior marginal probability 
#' e.g. $P(X(t)|\pmb{E})$ forall $X(t)\in\Pi'(t)$. This function uses oneslicesumprod and 
#' Pi-nodeorder/cptaddch/cptaddcht data structures to compute marginals for nodes in layer
#' $\pi_d(t)\in\Pi(t)$, then simply marginalise out other nodes in layer to get node marginal.
#' @param t = current time slice
#' @param Jb joint struct to initialise backwards sweep with (NULL if nothing to initialise), 
#' $p numeric vec of p, $n int vec of nodes n in $J$ in
#' order in which $J$ is composed, $t int vec of node t-slices 1..tend, $s list of state-labels (as 
#' indices) foreach node in $n$
#' @param indb - vector of indices to Jb$n/t for nodes to marginalise
#' @param Jf - like Jb but to "initialise" (elim at end) for forwards sweep due to 
#' backwards elimination order
#' @param indf - vector of indices to Jf$n/t for nodes to marginalise
#' @param Jout joint struct J[[t-slice,node n]] to write posterior marginals to same struct as Jinit
#' @param Jcpt[[1..tend,target node]] is the original joint populated by $cpt (same struct as J) 
#' @param Pinodeorder = is equivalent to lknodeorder wrapped by layer d, Pinodeorder[[d]]= defines 
#' the order of elimination of nodes in this layer given results from prev layers; $n = node id, 
#' $k = dt-slice
#' @param Picptaddch[[layer d]]  contains list of vec of nodeinds corresponding to CPT's to add for
#' multiply foreach elimination step in Pibknodeorder[[d]]
#' @param Picptaddcht[[layer d]]- ditto but for node dt-slice
#' @return Jout joint structure with node marginals forall nodes at slice t updated
#' @export
#' @examples none
#' @name oneslicemargs
oneslicemargs = function(t,Jb,indb,Jf,indf,Jout,Jcpt,Pinodeorder,Picptaddch,Picptaddcht) {
  # Iterate through all layers d=1 to D (R indexing starts at 1)
  for(i in 1:length(Pinodeorder)) {
    cat('Pinodeorder iteration ',i,'\n')
    J = oneslicesumprod(t,Jb=Jb,indb=indb,Jf=Jf,indf=indf,Jcpt,Pinodeorder[[i]],
                        Picptaddch[[i]],Picptaddcht[[i]])
    cat('J$n=',J$n,'J$t=',J$t,'\n')
    # Get marginals foreach node in layer d
    nt = data.frame(n=J$n,t=J$t)
    if(nrow(nt)>1) {
      for(j in 1:nrow(nt)) {
        Jout[[nt$t[j],nt$n[j]]] = marginalise(J$p,J$n,J$t,J$s,nt$n[-j],nt$t[-j])
        # cat('t=',nt$t[j],' n=',nt$n[j],'p=',Jout[[nt$t[j],nt$n[j]]]$p,'\n')
      }
    } else {
      Jout[[J$t,J$n]] = J
      # print(J)
    }
    Jout[[t,i]]$p = Jout[[t,i]]$p / sum(Jout[[t,i]]$p) # Normalise = P(\theta\cap Y)/P(Y)
  }
  return(Jout)
}
# # ---------------------------
# t1 = Sys.time()
# Jf = vector('list',length=tend) # Stores Jf to use at time slice t, which is actually $\Pi(t-1)$
# Jb = vector('list',length=tend)
# indf = NULL
# for(t in 1:5) { 
#   cat('t=',t,'-----------------------------------------------------------------\n')
#   cat('Jf has n=',Jf[[t]]$n,' t=',Jf[[t]]$t,'\n')
#   indf = which(Jf[[t]]$t != t)
#   J = oneslicemargs(t=t,Jb=NULL,indb=NULL,Jf=Jf[[t]],indf=indf,Jout=J,Jcpt,
#                      Pinodeorder,Picptaddch,Picptaddcht) # Compute marginals forall X(t)
#   if(t!=tend) {
#     Jf[[t+1]] = oneslicesumprod(t,Jb=NULL,indb=NULL,Jf=Jf[[t]],indf=indf,Jcpt,
#                                 lknodeorder,lkcptaddch,lkcptaddcht) #Compute P(\pi(t))
#   #   Jcpt = updJcpt(nnames,statenames,t=t,Jb=Jb,Jf=Jf,J,
#   #                  Jcpt,lknodeorder=lknodeorder,lkcptaddch,lkcptaddcht) # Update CPT at t+1
#   }
# }
# cat('Took ',Sys.time()-t1,'s\n')

# ----------------------------
#' @title oneslicemargs2
#' @description For a specified slice e.g. $\Pi'(t)$, find the posterior marginal probability 
#' e.g. $P(X(t)|\pmb{E})$ forall $X(t)\in\Pi'(t)$. This function is NAIVE, uses oneslicesumprod and 
#' nodeorder/cptaddch/cptaddcht data structures to compute marginals for one node at a time in 
#' $\Pi'(t)$. 
#' @param t = current time slice
#' @param Jb joint struct to initialise backwards sweep with (NULL if nothing to initialise), 
#' $p numeric vec of p, $n int vec of nodes n in $J$ in
#' order in which $J$ is composed, $t int vec of node t-slices 1..tend, $s list of state-labels (as 
#' indices) foreach node in $n$
#' @param indb - vector of indices to Jb$n/t for nodes to marginalise
#' @param Jf - like Jb but to "initialise" (elim at end) for forwards sweep due to 
#' backwards elimination order
#' @param indf - vector of indices to Jf$n/t for nodes to marginalise
#' @param Jout joint struct J[[t-slice,node n]] to write posterior marginals to same struct as Jinit
#' @param Jcpt[[1..tend,target node]] is the original joint populated by $cpt (same struct as J) 
#' @param nodeorder = is equivalent to lknodeorder wrapped by layer d, nodeorder[[d]]= defines 
#' the order of elimination of nodes in this layer given results from prev layers; $n = node id, 
#' $k = dt-slice
#' @param cptaddch[[layer d]]  contains list of vec of nodeinds corresponding to CPT's to add for
#' multiply foreach elimination step in nodeorder[[d]]
#' @param cptaddcht[[layer d]]- ditto but for node dt-slice
#' @return Jout joint structure with node marginals forall nodes at slice t updated
#' @export
#' @examples none
#' @name oneslicemargs2
oneslicemargs2 = function(t,Jb,indb,Jf,indf,Jout,Jcpt,nodeorder,cptaddch,cptaddcht) {
  # dim(J)
  for(i in 1:length(cptaddch)) {
    Jout[[t,i]] = oneslicesumprod(t,Jb=Jb,indb=indb,Jf=Jf,indf=indf,Jcpt,nodeorder[[i]],
                        cptaddch[[i]],cptaddcht[[i]])
    # print(Jout[[t,i]])
    Jout[[t,i]]$p = Jout[[t,i]]$p / sum(Jout[[t,i]]$p) # Normalise = P(\theta\cap Y)/P(Y)
  }
  return(Jout)
}

# --------------------------------------------------------------------------------------------------
#' @title updJcpt
#' @description Non-homogeneous DBN custom function for the seagrass model to update Jcpt given
#' J outcomes at previous time slice. Dynamically updates Net change Shoot Density CPT in Jcpt using
#' the estimate of the net change shoot density at time $t+1$, $P(\hat{S})$, found using variable
#' elimination. Based off net gain in shoot density (going from lower states to higher states) 
#' $\delta=P(\hat{S_{t+1}})-P(S_t)$. 3 options for updating CPT: a) overwrite entire CPT, b) write
#' in zero-zero-zero column of CPT, c) write pgain only to zero-zero-zero column of cpt.
#' @param nnames = char vector of DBN node names
#' @param statenames (char vec odf node state names),
#' @param t = current time slice
#' @param Jb[[1..tend]] joint struct to initialise backwards sweep with. $p numeric vec of p, $n int
#' vec of nodes n in $J$ in order in which $J$ is composed, $t int vec of node t-slices 1..tend, $s
#'  list of state-labels (as indices) foreach node in $n$
#' @param Jf[[1..tend]] - like Jb but to "initialise" (elim at end) for forwards sweep due to 
#' backwards elimination order
#' @param Jt - joint structure; we use Jt[[t,n]] = computed ACTUAL node marginal for updating CPTs.
#' Jt[[t,n]] $p numeric vector of p, $n int vec of nodes n in $J$ in order in which $J$ is composed,
#' $t int vec of node t-slices 1..tend, $s list of state labels (as indices) foreach node in $n$.
#' @param Jcpt[[1..tend,target node]] is the original joint populated by $cpt (same struct as J) 
#' @param lknodeorder = $\Pi'(t)$ defines the order of elimination of nodes in the slice in order
#' to get slice marginals e.g. $P(\pi'(t)|\pmb{E})$. lknodeorder $n = node id, $k = dt-slice
#' @param lkcptaddch = list of vec of nodeinds corresponding to CPT's to add for multiply foreach
#' elimination step in lknodeorder
#' @param lkcptaddcht = same as lkcptaddch but for node dt-slice
#' @param W = exogeneous factors list structure, W$Seed[1:tend] = adjusted count of #months of ploss
#' @param evid_t[[1..tend,node]][vector of len numstates containing node marginals] - used only
#' to suppress updating CPT if evidence is supplied for that node
#' @return list containing updated CPT structure Jcpt and updated exogeneous factors W
#' @export
#' @examples none
#' @name updJcpt
updJcpt = function(nnames,statenames,t,Jb,Jf,Jt,Jcpt,lknodeorder,lkcptaddch,lkcptaddcht,W,evid_t){
  ncind = grep('Net_Change_',nnames) # Index to net_change_ nodes
  nrind = grep('Realised_',nnames) # Index to realised shoot density/biomass nodes
  nlind = grep('Accumulated_Light',nnames) # Index to accumulated_light node
  nseed = which('Seed_Density'==nnames) # Index to seed density node
  nphys = which('Physiological_Status_of_Plants'==nnames)
  Jo = NULL # J structure storing oneslicesumprod result prior to marginalisation for shoot den/seed
  # ------------------------------------------------------------------------------------------------
  # Update net change shoot density
  for(i in 1:length(ncind)) { # foreach net_change/realised_ combo
    J=NULL
    n = ncind[i] # index to net change node
    nr = nrind[i] # index to shoot density/biomass node
    # Step 1: Estimate $P(\hat{S(t+1)})$ (realised shoot density)
    # if(!is.null(Jb)) { # Nodes to elim in Jb before getting started; e.g. $P(\pi(t+1)|E)$
    #   J = marginalise(Jb$p,Jb$n,Jb$t,Jb$s,Jb$n[indb],Jb$t[indb])
    # }
    # J = oneslicesumprod(t+1,J,Jcpt,lknodeorder,lkcptaddch,lkcptaddcht) # Eliminate t+1
    # if(!is.null(Jf)) { # Eliminate any nodes in $\P(pi(t-1)|\pmb{E})$
    #   J = multiply(J$p,Jf$p,J$n,Jf$n,J$t,Jf$t,J$s,Jf$s)
    #   J = marginalise(J$p,J$n,J$t,J$s,Jf$n[indf],Jf$t[indf])
    # }
    Jo = oneslicesumprod(t+1,Jb=Jb[[t+1]],indb=which(Jb[[t+1]]$t!=t+1),Jf=Jf[[t+1]],
                        indf=which(Jf[[t+1]]$t!=t+1),Jcpt,lknodeorder,lkcptaddch,lkcptaddcht)
    ind = which(Jo$n==nr & Jo$t==t+1) #Marginalise out other nodes to get estimate $P(\hat{S(t+1)})$
    J = marginalise(Jo$p,Jo$n,Jo$t,Jo$s,Jo$n[-ind],Jo$t[-ind])
    # cat('step 1 done, J=\n')
    # cat('ShootD: step 1 done, J$p=',J$p,'n=',J$n,'t=',J$t,'s=',J$s[[1]],'\n')
    # print(J)
    # Step 2: Find Pgain using $P(\hat{S}(t+1))-P(S(t))$ (diff in realised shoot d between t, t+1)
    ns = length(statenames[[n]])
    xst = vector('numeric',ns) # initialise with zeros
    xs = xst
    xst[J$s[[1]]] = J$p # write p to xst and xs according to state indices in case CPT is filted
    xs[Jt[[t,nr]]$s[[1]]] = Jt[[t,nr]]$p
    pgain = sum(c(3,2,1,0) * (xst - xs)) /2 # ASSUMPTION halve the gain works well empirically
    W$shootpgain[t+1] = pgain
    pgain_ = pgain
    xst_ = xst
    # cat('step 2 done, pgain=',pgain,'\n')
    # Step 3: distribute pgain over states
    pgain2 = 0
    if(pgain>0 && length(evid_t[[t+1,n]])==0) {
      for(j in ns:1) { # for j in z,l,m,h
        xst[j] = xs[j] - min(pgain,xs[j]) + pgain2
        pgain2 = min(pgain,xs[j])
        pgain = pgain - pgain2
        if(j==1)
          xst[j] = min(xst[j] + pgain,1)
      }
      #normalise as small errs ~1e-14 accumulate to collapse p-->0 in ~48 t-slices
      xst = xst / sum(xst)
      # cat('step 3 done, xst=',xst,'\n')
      # Step 4: write to CPT (of net change shoot density)
      # Jcpt should have n=12,23,24 and we are interested in column corresponding to state combo 4,4,4
      # Option a: overwrite entire CPT
      Jcpt[[t+1,n]]$p = xst
      Jcpt[[t+1,n]]$n = n
      Jcpt[[t+1,n]]$t = t+1
      Jcpt[[t+1,n]]$s = list(1:ns)
      # cat('option a done\n')
      # print(Jcpt[[t+1,n]])
      # Option b: overwrite zero-zero-zero column of CPT corresponding to stateid=4 for net chg node
      # Assumes order of nodes remains constant (since CPT for that node)
      # nsvect = sapply(Jcpt[[t+1,n]]$s,function(X)length(X))
      # if(max(Jcpt[[t+1,n]]$s[[length(nsvect)]])==4 && nsvect[length(nsvect)]==4) {
      # Only apply if zero state is still present stateid=4 for netchg AND 4 states for net chg node
      # Jcpt[[t+1,n]]$p[(prod(nsvect)-nsvect[length(nsvect)]+1):prod(nsvect)] = xst
      # cat('option b done\n')
      # print(Jcpt[[t+1,n]])
      # Option c: overwrite zero-zero-zero column with pgain only - same as b, just diff overwrite
      # Jcpt[[t+1,n]]$p[(prod(nsvect)-nsvect[length(nsvect)]+1):prod(nsvect)] = c(0,0,pgain_,1-pgain_)
      # cat('option c done\n')
      # print(Jcpt[[t+1,n]])
      # }
    } else {
      #return(Jcpt) # Make no change to CPT if pgain is negative - absorb to zero
    } # end if pgain>0
    
    if(xst_[ns]>=0.999){ # zero is last state, so index is ns=#states
      # Cut to seed p_zero=1 if p_Zero predicted at t+1 >=0.995
      Jcpt[[t+1,n]]$p = vector('numeric',ns)
      Jcpt[[t+1,n]]$p[ns] = 1
      Jcpt[[t+1,n]]$n = n
      Jcpt[[t+1,n]]$t = t+1
      Jcpt[[t+1,n]]$s = list(1:ns)
    }
  } # end foreach net change x/realised x combo
  
  # ------------------------------------------------------------------------------------------------
  # Update seed density if Halophila or Zostera
  if(any(Jcpt[[1,which(nnames %in% 'Genera_Presence')]]$s[[1]] %in% c(1,11)) && length(evid_t[[t+1,nseed]])==0) {
    # It is Halophila/Zostera, lets do this!
    J=NULL
    n = nseed # index to seed density node
    # Step 1: Estimate $P(\hat{S(t+1)})$ (seed density)
    # Use oneslicesumprod calculated for shoot density since that, seed and phys are all link nodes
    ind = which(Jo$n==n & Jo$t==t+1) #Marginalise out other nodes to get estimate $P(\hat{S(t+1)})$
    J = marginalise(Jo$p,Jo$n,Jo$t,Jo$s,Jo$n[-ind],Jo$t[-ind])
    # cat('step 1 done, J$p=',J$p,'n=',J$n,'t=',J$t,'s=',J$s[[1]],'\n')
    # print(Jt[[t,n]])
    # Step 2: Find Ploss using $P(\hat{S}(t+1))-P(S(t))$ (diff in seed d between t, t+1)
    ns = length(statenames[[n]])
    xst = vector('numeric',ns) # initialise with zeros
    xs = xst
    xst[J$s[[1]]] = J$p # write p to xst and xs according to state indices in case CPT is filted
    xst_ = xst # store (t+1) xs for cut-to-zero filtering below
    xs[Jt[[t,n]]$s[[1]]] = Jt[[t,n]]$p
    ploss = sum(c(0,1,2) * (xst - xs)) 
    W$seedploss[t+1] = ploss
    # cat('Pre-scale ploss = ',ploss,'\n')
    ploss = ploss * min(W$seed[t]+1,6)^2/6^2
    ploss_ = ploss
    W$seedplossadj[t+1] = ploss
    # cat('step 2 done, ploss=',ploss,'\n')
    # Step 3: distribute pgain over states
    ploss2 = 0
    if(ploss>1e-12) { 
      for(j in c(1,2,3)) { # for j in l,h,z
        xst[j] = xs[j] - min(ploss,xs[j]) + ploss2
        ploss2 = min(ploss,xs[j])
        ploss = ploss - ploss2
        if(j==3)
          xst[j] = min(xst[j] + ploss,1)
      }
      
      xst = xst/sum(xst)#normalise as small errs ~1e-14 accumulate to collapse p-->0 in ~48 t-slices
      # cat('step 3 done, xst=',xst,'\n')
      # Step 4: write to CPT (of net change shoot density)
      # Jcpt should have n=12,23,24; we are interested in column corresponding to state combo 4,4,4
      # Option a: overwrite entire CPT
      Jcpt[[t+1,n]]$p = xst
      Jcpt[[t+1,n]]$n = n
      Jcpt[[t+1,n]]$t = t+1
      Jcpt[[t+1,n]]$s = list(1:ns)
      W$seed[t+1] = min(6,W$seed[t]+1)
    } else {
      W$seed[t+1] = W$seed[t] # if no loss, and no reset (see below), set same as previous t-slice
    } 
    # if(abs(ploss_)<1e-2 && W$seed[t]>0){#Case of 'no change' but previously loss-so continue count
    #   W$seed[t+1] = min(12,W$seed[t]+1)
    # } else 
    if(ploss_<=-1e-12 && Jcpt[[t+1,nlind]]$s[[2]]!=2) { # Use 1e-12 due to floating point errors
      # (since ploss adjusted by W$seed); AND don't reset if still dredging at t+1
      # (light is still p=1 i.e. below_sat)
      # cat(Jcpt[[t+1,nlind]]$s[[2]],' aa\n')
      W$seed[t+1] = 0 #max(0,W$seed[t]-1) # Reset consecutive loss months to zero
    } # end if ploss>0
    
    if(xst_[ns]>=0.999){ # zero is last state, so index is ns=#states
      # Cut to seed p_zero=1 if p_Zero predicted at t+1 >=0.995 
      Jcpt[[t+1,n]]$p = vector('numeric',ns)
      Jcpt[[t+1,n]]$p[ns] = 1
      Jcpt[[t+1,n]]$n = n
      Jcpt[[t+1,n]]$t = t+1
      Jcpt[[t+1,n]]$s = list(1:ns)
    }
  }
  # ------------------------------------------------------------------------------------------------
  # Update phys status if Amphibolis
  if(any(Jcpt[[1,which(nnames %in% 'Genera_Presence')]]$s[[1]]%in%c(2)) && length(evid_t[[t+1,nphys]])==0) {
    # It is Amphibolis, lets do this!
    J=NULL
    n = nphys # index to seed density node
    # Step 1: Estimate $P(\hat{S(t+1)})$ (physiological status)
    # Use oneslicesumprod calculated for shoot density since that, seed and phys are all link nodes
    ind = which(Jo$n==n & Jo$t==t+1) #Marginalise out other nodes to get estimate $P(\hat{S(t+1)})$
    J = marginalise(Jo$p,Jo$n,Jo$t,Jo$s,Jo$n[-ind],Jo$t[-ind])
    # cat('step 1 done, J=\n')
    # print(J)
    # Step 2: Find Ploss using $P(\hat{S}(t+1))-P(S(t))$ (diff in phys status between t, t+1)
    ns = length(statenames[[n]])
    xst = vector('numeric',ns) # initialise with zeros
    xs = xst
    xst[J$s[[1]]] = J$p # write p to xst and xs according to state indices in case CPT is filted
    xs[Jt[[t,n]]$s[[1]]] = Jt[[t,n]]$p
    ploss = sum(c(0,1,2) * (xst - xs)) 
    W$physploss[t+1] = ploss
    ploss = ploss * min(W$phys[t]+1,2)^2/4
    W$physplossadj[t+1] = ploss
    ploss_ = ploss
    # cat('step 2 done, ploss=',ploss,'\n')
    # Step 3: distribute pgain over states
    ploss2 = 0
    if(ploss>0) { 
      for(j in c(1,2,3)) { # for j in l,h,z
        xst[j] = xs[j] - min(ploss,xs[j]) + ploss2
        ploss2 = min(ploss,xs[j])
        ploss = ploss - ploss2
        if(j==3)
          xst[j] = min(xst[j] + ploss,1)
      }
      
      xst = xst/sum(xst)#normalise as small errs ~1e-14 accumulate to collapse p-->0 in ~48 t-slices
      # cat('step 3 done, xst=',xst,'\n')
      # Step 4: write to CPT (of net change shoot density)
      # Jcpt should have n=12,23,24; we are interested in column corresponding to state combo 4,4,4
      # Option a: overwrite entire CPT
      Jcpt[[t+1,n]]$p = xst
      Jcpt[[t+1,n]]$n = n
      Jcpt[[t+1,n]]$t = t+1
      Jcpt[[t+1,n]]$s = list(1:ns)
      W$phys[t+1] = min(2,W$phys[t]+1)
    } 
    # if(abs(ploss_)<1e-2 && W$phys[t]>0){#Case of 'no change' but previously loss-so continue count
    #   W$phys[t+1] = min(12,W$phys[t]+1)
    # } else 
    if(ploss_<=1e-12) { # Use 1e-12 due to floating point errors (since ploss adjusted by W$phys)
      W$phys[t+1] = 0 #max(0,W$phys[t]-1) # Reset consecutive loss months to zero
    } # end if ploss>0
  }
  return(list(Jcpt=Jcpt,W=W))
}

# --------------------------------------------------------------------------------------------------
#' @title fwdonly
#' @description Non-homogeneous DBN forwards only inference. Uses oneslicemargs2, oneslicesumprod,
#' and updJcpt to calculate marginals for X(t) forall t, compute slice partial results for $\pi(t)$,
#' and update non-homogeneous CPTs respectively. Compatible with naive oneslicemargs as well as 
#' d-layer based oneslicemargs.
#' @param tend = index to last t-slice in simulation window 1:tend
#' @param Jcpt[[1..tend,target node]] is the original joint populated by $cpt (same struct as J) 
#' @param J = J[[1..tend,target node]] - Joint structure for inference with $p numeric vector of p,
#' $n int vec of nodes n in $J$ in order in which $J$ is composed, $t int vec of node t-slices
#' 1..tend, $s list of state labels (as indices) foreach node in $n$.
#' @param Pinodeorder = is equivalent to lknodeorder wrapped for a given layer d OR for a given
#' node. Pinodeorder[[d]]= defines the order of elimination of nodes in this layer given results 
#' from prev layers; $n = node id, $k = dt-slice
#' @param Picptaddch[[layer d]]  contains list of vec of nodeinds corresponding to CPT's to add for
#' multiply foreach elimination step in Pibknodeorder[[d]]
#' @param Picptaddcht[[layer d]]- ditto but for node dt-slice
#' @param lknodeorder = $\Pi'(t)$ defines the order of elimination of nodes in the slice in order
#' to get slice marginals e.g. $P(\pi'(t)|\pmb{E})$. lknodeorder $n = node id, $k = dt-slice
#' @param lkcptaddch = list of vec of nodeinds corresponding to CPT's to add for multiply foreach
#' elimination step in lknodeorder
#' @param lkcptaddcht = same as lkcptaddch but for node dt-slice
#' @param nnames = char vector of DBN node names
#' @param statenames (char vec odf node state names),
#' @param W = exogeneous factors list structure, W$Seed[1:tend] = adjusted count of #months of ploss
#' @param evid_t[[1..tend,node]][vector of len numstates containing node marginals] - used only
#' to suppress updating CPT updJcpt if evidence is supplied for that node
#' @return list containing updated $J[[t,node]] of node marginals, $Jcpt[[t,node]] node CPTs,
#' $Jf[[1..tend]] $\pi(t)$ partial results to be used at time t for forwards sweep, 
#' $Jb[[1..tend]] ditto for backwards sweep, W$seed[1:tend] exogeneous factors, seed dens losscount.
#' @export
#' @examples none
#' @name fwdonly
# Note does not work for $\Pi'(t)$, only for $\Pi(t)$
fwdonly = function(tend,Jcpt,J,Pinodeorder,Picptaddch,Picptaddcht,
                   lknodeorder,lkcptaddch,lkcptaddcht,nnames,statenames,W,evid_t) {
  t1 = Sys.time()
  Jf = vector('list',length=tend) # Stores Jf to use at time slice t, which is actually $\Pi(t-1)$
  Jb = vector('list',length=tend)
  indf = NULL
  for(t in 1:tend) {
    cat('t=',t,'-----------------------------------------------------------------\n')
    indf = which(Jf[[t]]$t != t)
    J = oneslicemargs2(t=t,Jb=NULL,indb=NULL,Jf=Jf[[t]],indf=indf,Jout=J,Jcpt,
                      Pinodeorder,Picptaddch,Picptaddcht) # Compute marginals forall X(t)
    if(t!=tend) {
      Jf[[t+1]] = oneslicesumprod(t,Jb=NULL,indb=NULL,Jf=Jf[[t]],indf=indf,Jcpt,
                                lknodeorder,lkcptaddch,lkcptaddcht) #Compute P(\pi(t))
      X = updJcpt(nnames,statenames,t=t,Jb=Jb,Jf=Jf,J,
                     Jcpt,lknodeorder=lknodeorder,lkcptaddch,lkcptaddcht,W,evid_t)#Update CPT at t+1
      Jcpt = X$Jcpt
      W = X$W
    }
  }
  cat('Took ',Sys.time()-t1,'s\n')
  return(list(J=J,Jcpt=Jcpt,Jf=Jf,Jb=Jb,W=W))
}


# --------------------------------------------------------------------------------------------------
#' @title fwdback
#' @description Non-homogeneous DBN forwards and backwards inference. Uses oneslicemargs2, 
#' oneslicesumprod, and updJcpt to calculate marginals for X(t) forall t, compute slice partial 
#' results for $\pi(t)$, and update non-homogeneous CPTs respectively. Compatible with naive 
#' oneslicemargs as well as d-layer based oneslicemargs. Assumes default model for t=1 and for 
#' t=2:tend has already been put into Jcpt. Also, updates CPTs using updJcpt.
#' @param tend = index to last t-slice in simulation window 1:tend
#' @param Jcpt[[1..tend,target node]] is the original joint populated by $cpt (same struct as J) 
#' @param J = J[[1..tend,target node]] - Joint structure for inference with $p numeric vector of p,
#' $n int vec of nodes n in $J$ in order in which $J$ is composed, $t int vec of node t-slices
#' 1..tend, $s list of state labels (as indices) foreach node in $n$.
#' @param Pinodeorder = is equivalent to lknodeorder wrapped for a given layer d OR for a given
#' node. Pinodeorder[[d]]= defines the order of elimination of nodes in this layer given results 
#' from prev layers; $n = node id, $k = dt-slice
#' @param Picptaddch[[layer d]]  contains list of vec of nodeinds corresponding to CPT's to add for
#' multiply foreach elimination step in Pibknodeorder[[d]]
#' @param Picptaddcht[[layer d]]- ditto but for node dt-slice
#' @param lknodeorder = $\Pi'(t)$ defines the order of elimination of nodes in the slice in order
#' to get slice marginals e.g. $P(\pi'(t)|\pmb{E})$. lknodeorder $n = node id, $k = dt-slice
#' @param lkcptaddch = list of vec of nodeinds corresponding to CPT's to add for multiply foreach
#' elimination step in lknodeorder
#' @param lkcptaddcht = same as lkcptaddch but for node dt-slice
#' @param nnames = char vector of DBN node names
#' @param statenames (char vec odf node state names),
#' @param W = exogeneous factors list structure, W$Seed[1:tend] = adjusted count of #months of ploss
#' @param evid_t[[1..tend,node]][vector of len numstates containing node marginals] - used only
#' to suppress updating CPT updJcpt if evidence is supplied for that node
#' @return list containing updated $J[[t,node]] of node marginals, $Jcpt[[t,node]] node CPTs,
#' $Jf[[1..tend]] $\pi(t)$ partial results to be used at time t for forwards sweep, 
#' $Jb[[1..tend]] ditto for backwards sweep, W$seed[1:tend] exogeneous factors, seed dens losscount.
#' @export
#' @examples none
#' @name fwdback
# Note does not work for $\Pi'(t)$, only for $\Pi(t)$
fwdback = function(tend,Jcpt,J,Pinodeorder,Picptaddch,Picptaddcht,
                   lknodeorder,lkcptaddch,lkcptaddcht,nnames,statenames,W,evid_t) {
  t1 = Sys.time()
  Jf = vector('list',length=tend) # Stores Jf to use at time slice t, which is actually $\Pi(t-1)$
  Jb = vector('list',length=tend)
  indf = NULL
  indb = NULL
  
  # Pre-compute
  # for(t in 1:(tend-1)) { # Forwards sweep pre-compute
  #   indf = which(Jf[[t]]$t != t)
  #   Jf[[t+1]] = oneslicesumprod(t,Jb=NULL,indb=NULL,Jf=Jf[[t]],indf=indf,Jcpt,
  #                               lknodeorder,lkcptaddch,lkcptaddcht) #Compute partial fwd P_f(\pi(t))
  # }
  for(t in tend:2) { # Backwards sweep pre-compute partial accumulated backwards P_b(\pi(tend:t))
    indb = which(Jb[[t]]$t != t)
    Jb[[t-1]] = oneslicesumprod(t,Jb=Jb[[t]],indb=indb,Jf=NULL,indf=NULL,Jcpt,
                                lknodeorder,lkcptaddch,lkcptaddcht) #Compute partial bwd P_b(\pi(t))
  }
  # Iterate
  for(t in 1:tend) {
    cat('t=',t,'-----------------------------------------------------------------\n')
    indf = which(Jf[[t]]$t != t)
    indb = which(Jb[[t]]$t != t)
    J = oneslicemargs2(t=t,Jb=Jb[[t]],indb=indb,Jf=Jf[[t]],indf=indf,Jout=J,Jcpt,
                       Pinodeorder,Picptaddch,Picptaddcht) # Compute marginals forall X(t)
    if(t!=tend) {
      Jf[[t+1]] = oneslicesumprod(t,Jb=NULL,indb=NULL,Jf=Jf[[t]],indf=indf,Jcpt,
                                  lknodeorder,lkcptaddch,lkcptaddcht) #Compute P(\pi(1:t))
      X = updJcpt(nnames,statenames,t=t,Jb=Jb,Jf=Jf,J,
                     Jcpt,lknodeorder=lknodeorder,lkcptaddch,lkcptaddcht,W,evid_t)#Update CPT at t+1
      Jcpt = X$Jcpt
      W = X$W
    }
  }
  cat('Took ',Sys.time()-t1,'s\n')
  return(list(J=J,Jcpt=Jcpt,Jf=Jf,Jb=Jb,W=W))
}


# --------------------------------------------------------------------------------------------------
#' @title initdbfivefdsPi
#' @description Initialise DBFIVE From DbnSmiler object and evidence evid for Pi based inferencing 
#' @param dbn - DBNSMILER object 
#' @param E - DBNSMILER evidence structure (dbn$evid)
#' @param maxdt - max dt-slice (=2 for Markov system)
#' @return list containing:
#' \itemize{
#'  \item $nnames (char vec of node names),
#'  \item $statenames (char vec odf node state names),
#'  \item $parnodeind ([[1..maxdt,node]][parent node ind to nnames]),
#'  \item $parnodetslice ([[1..maxdt,node]][par k for $\Pi(t-k)$]), 
#'  \item $childnodeind ([[1..maxdt,node]][child node ind to nnames])
#'  \item $childnodetslice ([[1..maxdt,node]][child k for $\Pi(t+k)$]), 
#'  \item $numstates (vec of number of node states),
#'  \item $cpt ([[1..maxdt,node]][numeric vector of CPT])
#'  \item $maxdt = (int) max# dt-slices (e.g. Markov network, maxdt = 1; BN, maxdt=0)
#'  \item $tend = index to last t-slice in simulation window 1:tend
#'  \item $evid_t [[1..tend,node]][vector of len numstates containing node marginals]
#'  \item $nodecost[1..maxdt,node] containing cost value $OC=D+C$
#'  \item $nodedist[1..maxdt,node] containing distance $D$
#'  \item $parpars
#'    \itemize{
#'      \item $n[[k=1..maxdt for $\Pi(t-k)$,node]] vector of node inds
#'      \item $t[[k=1..maxdt for $\Pi(t-k)$,node]] vector of node dt-slice t-k
#'    }
#'  \item $nodeorder 
#'    \itemize{
#'      \item $n[[k=1..maxdt for $\Pi(t-k)$,target node]] vector of node inds in order to eliminate
#'      \item $t[[k=1..maxdt for $\Pi(t-k)$,target node]] vector of node dt-slice t-k
#'    }
#'  \item $nodeorder2[[foreach node]] (which is nodeorder[[maxdt,]] reformatted in nodeorder2[[n]]$n
#'  and nodeorder2[[n]]$k storing vectors of node ids and dt-slice t-k respectively for fwdonly()
#'  \item $linknodes = vector of node ids for link nodes $\pi(t)$
#'  \item $lknodeorder - like nodeorder but for link nodes $\pi(t)$
#'    \itemize{
#'      \item $n[[k=1..maxdt for $\Pi(t-k)$,target node]] vector of node inds in order to eliminate
#'      \item $t[[k=1..maxdt for $\Pi(t-k)$,target node]] vector of node dt-slice t-k
#'    }
#'  \item $nodeordch[[k=1..maxdt for $\Pi(t-k)$,target node]][[foreach nodeorder node]][vec of 
#'  childnodeind]
#'  \item $nodeordcht[[k=1..maxdt for $\Pi(t-k)$,target node]][[foreach nodeorder node]][vec of 
#'  child dt-slice k as t+k]
#'  \item $cptaddch[[k=1..maxdt for $\Pi(t-k)$,target node]][[foreach nodeorder node]][vec of 
#'  nodeind corresponding to CPT's to add for multiply step]
#'  \item $cptaddcht[[k=1..maxdt for $\Pi(t-k)$,target node]][[foreach nodeorder node]][vec of node
#' dt-slice corresponding to CPT's to add for multiply step]
#'  \item $lkcptaddch[[k=1..maxdt for $\Pi(t-k)$,target node]][[foreach nodeorder node]][vec of 
#'  nodeind corresponding to CPT's to add for multiply step] - for link nodes $\pi(t)$
#'  \item $lkcptaddcht[[k=1..maxdt for $\Pi(t-k)$,target node]][[foreach nodeorder node]][vec of 
#'  node dt-slice corresponding to CPT's to add for multiply step] - for link nodes $\pi(t)$
#' \item J[[1..tend,target node]] - Joint structure
#'    \itemize{
#'      \item $p = numeric vector of p
#'      \item $n = int vector of nodes n in $J$ in order in which $J$ is composed
#'      \item $t = int vector of node t-slices 1..tend
#'      \item $s = list of state-labels (as indices) foreach node in $n$
#'    }
#'  \item Jlink[[1..tend,target node]] is the joint for $\pi(t)$ (same struct as J)
#'  \item Jcpt[[1..tend,target node]] is the original joint populated by $cpt (same struct as J)
#'  \item nodeorder3 - used to find fbnodeorder etc., not used otherwise
#'  \item fbnodeorder[[foreach node]]$n and $k - same as nodeorder2 but for fwdback() inference
#'  \item fbcptaddch[[1,foreach node]] - like cptaddch but for fwdback() inference 
#'  \item fbcptaddcht[[1,foreach node]] - like cptaddcht but for fwdback() inference
#' } 
#' @export
#' @examples none
#' @name initdbfivefdsPi
initdbfivefdsPi = function(dbn,E,tend,maxdt) {
  # ----------------------- Conversion from DBNSMILER to DFIVE -----------------------
  x = dbn2DFnodes(dbn,maxdt) # Convert DBN structure to DFIVE structure # UP TO HERE
  nnames = x$nnames
  statenames = x$statenames
  parnodeind = x$parnodeind
  parnodetslice = x$parnodetslice
  childnodeind = x$childnodeind
  childnodetslice = x$childnodetslice
  numstates = x$numstates
  cpt = x$cpt
  
  evid_t = dbn2DFevid(dbn,E,nnames) # Convert DBN evid structure to DFIVE evid structure
  
  # ----------------------- Initialise DBFIVE data structures -----------------------
  x = calcnodecost(parnodeind,parnodetslice,childnodeind,childnodetslice,numstates,maxdt)
  nodecost = x$nodecost # Get nodecost
  nodedist = x$nodedist
  
  parpars = getparpars(parnodeind,parnodetslice,maxdt,maxdt) #Get parpars structure
  # Get nodeorder
  nodeorder = getnodeorder(parpars,childnodeind,childnodetslice,nodecost[1,],maxdt,nnames,NULL)
  nodeorder2 = vector('list',length(nodeorder$n[1,]))
  for(i in 1:length(nodeorder2)) {
    nodeorder2[[i]] = list(n = nodeorder$n[[1,i]], k=nodeorder$k[[1,i]])
  }
  
  # Filter childnodes to create list structures: nodeordch, nodeordcht and cptaddch/cptaddcht
  targetnodes = lapply(rep(1:length(nnames),times=2),function(X)X)
  targetnodes = matrix(targetnodes,nrow=maxdt,ncol=length(nnames),byrow=TRUE)
  x = filterchildnode(childnodeind,childnodetslice,nodeorder,targetnodes,maxdt)
  nodeordch = x$nodeordch # [[k,target node]][[foreach nodeorder node]][vec of child nodeind]
  nodeordcht = x$nodeordcht # ditto but child node dt-slice
  cptaddch = x$cptaddch # ditto but child nodeind cpt's to add for multply step
  cptaddcht = x$cptaddcht # ditto but cpt dt-slice
  
  # ----------------------- Initialise inference data structures -----------------------
  # Initialise joint structure J[[1..tend,target node]]
  J = vector('list',tend*length(nnames))
  J = lapply(J,function(X) list(p=numeric(0),n=integer(0),t=integer(0),s=integer(0)))
  dim(J) = c(tend,length(nnames))
  # Initialise and pre-load evidence, the 'filtered cpt' into J
  Jcpt = filtercpt(J,nnames,maxdt,parnodeind,parnodetslice,numstates,childnodeind,
                   childnodetslice,evid_t,cpt,tend)
  
  # Initialise linknodes, lknodeorder, lkctpaddch, lkcptaddcht,Pinodeorder, Picptaddch, Picptaddcht
  linknodes = which(sapply(childnodetslice[maxdt,],function(X)any(X>0))) # nodes with child at t+k
  nt = unlist(sapply(linknodes,function(X,nodeorder)nodeorder$n[[maxdt,X]],nodeorder=nodeorder)) + 
    IDTSCALE*unlist(sapply(linknodes,function(X,nodeorder)nodeorder$k[[maxdt,X]],nodeorder=nodeorder))
  nt = nt[!duplicated(nt) & !(nt %in% linknodes)] # not duplicated and not a link node
  nt = nt[order(nodecost[maxdt,floor(nt)])]
  lknodeorder=list(n=matrix(list(floor(nt)),1,1),
                   k=matrix(list(1/IDTSCALE*(nt-floor(nt))),1,1) )#lknodeorder equiv to nodeorder
  x = filterchildnode(childnodeind,childnodetslice,lknodeorder,matrix(list(linknodes),1,1),maxdt)
  lkcptaddch = x$cptaddch[[1,1]]
  lkcptaddcht = x$cptaddcht[[1,1]]
  lk = lknodeorder
  lknodeorder = list()
  lknodeorder$n = lk$n[[1]]
  lknodeorder$k = lk$k[[1]]
  
  # Setup fbnodeorder, fbcptaddch, fbcptaddcht for forwards-backwards inference
  nodeorder3 = list(n=list(),k=list())
  nodeorder3$n = lapply(1:ncol(nodeorder$n),function(X) setdiff(order(nodecost[maxdt,]),X) )
  nodeorder3$k = lapply(nodeorder3$n, function(X) X*0)
  dim(nodeorder3$n) = c(1,length(nodeorder3$n))
  dim(nodeorder3$k) = c(1,length(nodeorder3$k))
  x = filterchildnode(childnodeind,childnodetslice,nodeorder3,targetnodes,maxdt)
  fbcptaddch = x$cptaddch # ditto but child nodeind cpt's to add for multply step
  fbcptaddcht = x$cptaddcht # ditto but cpt dt-slice
  fbnodeorder = list()
  for(i in 1:ncol(nodeorder3$n)) {
    fbnodeorder[[i]] = list()
    fbnodeorder[[i]]$n = nodeorder3$n[[1,i]]
    fbnodeorder[[i]]$k = nodeorder3$k[[1,i]]
  }
  
  # Setup exogeneous factors to use when updJcpt - W
  W = list()
  W$seed = vector('numeric',tend) # variable storing #consecutive seed density loss t-slices
  W$seedploss =vector('numeric',tend) 
  W$seedplossadj = vector('numeric',tend) 
  W$phys = vector('numeric',tend) # variable storing #consecutive phys status loss t-slices
  W$physploss = vector('numeric',tend)
  W$physplossadj = vector('numeric',tend)
  W$shootpgain = vector('numeric',tend)
  
  return(
    list(nnames=nnames,statenames=statenames,parnodeind=parnodeind,parnodetslice=parnodetslice,
         childnodeind=childnodeind,childnodetslice=childnodetslice,numstates=numstates, cpt=cpt,
         evid_t=evid_t, nodecost=nodecost, nodedist=nodedist,parpars=parpars, 
         nodeorder=nodeorder, nodeorder2=nodeorder2, nodeordch=nodeordch,
         nodeordcht=nodeordcht, cptaddch=cptaddch, cptaddcht=cptaddcht, 
         J=J, Jcpt=Jcpt,W=W,
         linknodes=linknodes,lknodeorder=lknodeorder,lkcptaddch=lkcptaddch,lkcptaddcht=lkcptaddcht,
         nodeorder3=nodeorder3,fbnodeorder=fbnodeorder,fbcptaddch=fbcptaddch,fbcptaddcht=fbcptaddcht
         ))
}

#' @title initdbfiveevid
#' @description Initialise DBFIVE evid_t object, re-initialising J and Jcpt (all previous values
#' will be lost). 
#' @param evid_t[[t,node]][states] structure
#' @param I - DBFIVE data structure (for fields, see initdbfivefdsPi)
#' @param tend = index to last t-slice in simulation window 1:tend
#' @param maxdt - max dt-slice (=2 for Markov system)
#' @return Updated DBFIVE data structure I with J and Jcpt updated
#' @export
#' @examples none
#' @name initdbfivefdsPi (for fields, see initdbfivefdsPi)
initdbfiveevid = function(evid_t,I,tend,maxdt) {
  # Initialise joint structure J[[1..tend,target node]]
  J = vector('list',tend*length(I$nnames))
  J = lapply(J,function(X) list(p=numeric(0),n=integer(0),t=integer(0),s=integer(0)))
  dim(J) = c(tend,length(I$nnames))
  # Initialise and pre-load evidence, the 'filtered cpt' into J
  I$Jcpt = filtercpt(J,I$nnames,maxdt,I$parnodeind,I$parnodetslice,I$numstates,I$childnodeind,
                   I$childnodetslice,evid_t,I$cpt,tend)
  I$J = J
  I$evid_t = evid_t
  I$W$seed = 0*I$W$seed # initialise exogeneous factors W$seed to zero
  I$W$seedploss = 0*I$W$seedploss
  I$W$seedplossadj = 0*I$W$seedplossadj
  I$W$phys = 0*I$W$phys # initialise exogeneous factors W$phys to zero
  I$W$physploss = 0*I$W$physploss
  I$W$physplossadj = 0*I$W$physplossadj
  I$W$shootpgain = 0*I$W$shootpgain
  return(I)
}

# --------------------------------------------------------------------------------------------------
#' @title dbfiveJ2dbnsP
#' @description Convert marginal probabilities in DBFIVE J structure into DBNSMILER posterior P 
#' data frame for visualisation. 
#' @param dbn - DBNSMILER object 
#' @param scn - char scenario id
#' @param J[[1..tend,target node]] - Joint structure storing marginal posteriors, $p=numeric vector
#' of p, $n=int vec of nodes n in $J$ in order in which $J$ is composed, $t int vec of node t-slices
#' 1..tend, and $s=list of state-labels (as indices) foreach node in $n
#' @param nnames (char vec of node names)
#' @param statenames (char vec odf node state names)
#' @return DBNSMILER posterior P (data frame: node=char,state=char,t=char,p=numeric,scenario=int,
#'  pctile=numeric,xbar=numeric expected value for given percentile,ybar=observed numeric exp. val
#'  for given percentile, py = numeric observed probability)
#' @export
#' @examples none
#' @name dbfiveJ2dbnsP
dbfiveJ2dbnsP = function(dbn,scn,J,nnames,statenames) {
  # Convert J into dbn posterior structure for visualisation
  P = dbn$posterior
  P = P[-(1:nrow(P)),] # Clear current posterior
  ind = 1 # row index to P
  tvect = 1:tend + dbn$globdater$num[1] - 1 # Get datetime data from dbn structure
  tvect = num2chardate(tvect,dbn$datetype)
  for(j in 1:ncol(J)) { # foreach node
    for(i in 1:nrow(J)) { # foreach t-slice
      # s = c(J[[i,j]]$s[[1]], setdiff(1:length(statenames[[j]]), J[[i,j]]$s[[1]]) ) # state indices
      # sn = statenames[[j]][order(s)]
      # p = c(J[[i,j]]$p, rep(0.0,length(s)-length(J[[i,j]]$p)))[order(s)]
      p = rep(0.0,length(statenames[[j]]))
      p[J[[i,j]]$s[[1]]] = J[[i,j]]$p
      sn = statenames[[j]]
      
      for(k in 1:length(sn)) {
        P[ind,] = list(nnames[j],sn[k],tvect[i],p[k],scn,NA,NA,NA,NA)
        ind = ind + 1
      } # end foreach state
    } # end foreach t-slice
  } # end foreach node
  P = fastwmean(P,dbn,c(0,0.25,0.5,0.75,1)) # Get weighted mean
  return(P)
}


####################################################################################################
# >>>>>>>>>>>>>> ARCHIVE of
# --------------------------------------------------------------------------------------------------
#' @title checkmem DEPRECATED
#' @description simulate backwards inference (going from tend to t=1) by iterating over
#' maxdt slices to find the max memory required to store joint distribution J across all nodes.
#' Note this does not take into account $\lambda$ link nodes in runlkinf - only for runinf.
#' @param  nnames,maxdt,nodeorder,cptaddch,cptaddcht,parnodeind,parnodetslice
#' @return jmem[[k,node]][vector of #elements in J foreach nodeorder elim step]
#' @export
#' @examples none
#' @name checkmem
checkmem = function(nnames,maxdt,nodeorder,cptaddch,cptaddcht,parnodeind,parnodetslice) {
  jmem = vector('list',length(nnames)*maxdt)
  dim(jmem) = c(maxdt,length(nnames))
  for(j in 1:length(nnames)) { # foreach node
    jn = vector('integer',0)
    js = jn
    jnt = vector('numeric',0)
    for(i in 1:maxdt) { # foreach t-k scenario
      jmem[[i,j]] = vector('numeric',length(nodeorder$n[[i,j]]))
      for(l in 1:length(cptaddch[[i,j]])) { 
        if(length(cptaddch[[i,j]])>0) {
          n = cptaddch[[i,j]][[l]]#nodes to elim at step l ('L')
          np = unlist(sapply(parnodeind[i,n],function(X)X)) # par nodes of n
          nt = n + IDTSCALE*(maxdt-i+1-cptaddcht[[i,j]][[l]]) #t-slice(n)
        } 
        npt = np + IDTSCALE*(maxdt-i+1-unlist(sapply(parnodetslice[i,n],function(X)X)))#t-slice(np)
        jnt = c(jnt,nt,npt)
        jn = c(jn,n,np)
        jn = jn[!duplicated(jnt)] # Remove duplicates using jnt
        jnt = jnt[!duplicated(jnt)]
        js = numstates[jn]
        jmem[[i,j]][l] = prod(js) # stores 'worst case mem' before elimination
        # Eliminate node n
        mask = !(jnt %in% (nodeorder$n[[i,j]][l]+IDTSCALE*(maxdt-i+1-nodeorder$k[[i,j]][l])) )
        jn = jn[mask]
        jnt = jnt[mask]
      } # end foreach node to eliminate
    } # end foreach k
  } # end foreach node
  return(jmem)
}


# --------------------------------------------------------------------------------------------------
#' @title customPi - DEPRECATED (USER SHOULD SUPPLY OWN CUSTOMPI)
#' @description customPi - Create custom Pi' structure for Seagrass model. It does two things - 
#' (i) improve computational efficiency by assigning link nodes (see below) and (ii) enabling 
#' dynamic CPT updates for net change shoot density. <br/>
#' Create $\Pi'$ using heuristic $\Pi'$ = $\pi(t-1) \cup \Pi(t) \setminus \pi(t)$ for Markov DBNs
#' t-k where $k=1$. The rationale is that 'feedback' nodes $\pi(t-1)$ such as 'shoot density' have 
#' many child nodes in $t$; if these child nodes are eliminated first, that can greatly reduce the
#' computational burden in the multiply step.<br/>
#' @param nodeorder,ppn,ppt,maxdt,childnodeind,childnodetslice - see getnodeorder <br/>
#' @return nodeorder$n[[maxdt,numnodes]]=vector of nodeinds, $t[[maxdt,numnodes]]=vec of k+1
#' @export
#' @examples none
#' @name customPi
customPi = function(nodeorder,j,ppn,ppt,maxdt,childnodeind,childnodetslice,nnames) {
  K = as.integer(maxdt-1)
  
  # Create $\Pi'$ using heuristic $\Pi' = $\pi(t-1) \cup \Pi(t) \setminus \pi(t)$ for Markov t-1
  C=lapply(ppn,function(X,CNI,mdt) CNI[[mdt,X]],CNI=childnodeind,mdt=maxdt)#childnodeinds of ppn
  CT = mapply(function(X,Y) childnodetslice[[maxdt,X]]+Y,X=ppn,Y=ppt)#childnodetslices of ppn
  pi = vector('logical',length(ppn)) # mask - true if ppn/ppt belongs to \pi
  
  # $\Pi' = $\pi(t-1) \cup \Pi(t) \setminus \pi(t)$
  for(i in 1:length(ppn)) {
    pi[i] = any( sapply(1:length(CT[[i]]),function(X,ct,c,ppt,ppn,mask,i) 
      any(ct[X] != ppt[i] & ct[X]==ppt & c[X]==ppn),ct=CT[[i]],c=C[[i]],ppt=ppt,ppn=ppn,i=i) )
    if(is.na(pi[i]))
      pi[i] = FALSE # This happens when ppn[i] has no children, set as FALSE
  }
  
  # For k=1, get ppt==1 nodes, remove pi nodes and add pi nodes at t-1
  nodeorder$n[[2,j]] = c(ppn[ppt==1 & !pi],ppn[pi])
  nc = nodecost[nodeorder$n[[2,j]]]
  nodeorder$n[[2,j]] = nodeorder$n[[2,j]][order(nc)]
  nodeorder$k[[2,j]] = (K-c(ppt[ppt==1 & !pi],ppt[pi]-as.integer(1)))[order(nc)]#turn into k = K-ppt
  
  # For k=0, get ppt==2 nodes, drop first element = actual node - do not eliminate this!
  nodeorder$n[[1,j]] = ppn[ppt==2][-1]
  nc = nodecost[nodeorder$n[[1,j]]]
  nodeorder$n[[1,j]] = nodeorder$n[[1,j]][order(nc)]
  nodeorder$k[[1,j]] = (K+as.integer(1) - ppt[ppt==2][-1])[order(nc)]
  
  # Use special $\Pi'$ for dynamic CPT adaptation for net change target nodes so loss, recovery 
  # and net change and all child/grand-child of net change in same t-slice are put into the
  # previous "t-slice". Eg. CPT_p: p(zero,t)|p(zero,t-1) = f(loss(t-1),recov(t-1)). To do this,
  # wetake ppno/ppto and simply put all parents of net change, net change, and all children/grand-
  # children into t-1.
  if(substr(nnames[j],0,11)=='Net_Change_') {
    metricname = substr(nnames[j],12,nchar(nnames[j]))
    # Get chi-chi (child, grandchild etc. nodes) in the *same* t-slice - for net change, the nodes
    # are: Net Change, Loss in, Rate of Recovery in, and Realised
    inds = c(j, which(nnames %in% c(paste('Loss_in_',metricname,sep=''),paste('Rate_of_Recovery_in_',metricname,sep=''),
                                    paste('Realised_',metricname,sep=''))) )
    # For k=1, get ppt==1 nodes
    nodeorder$n[[2,j]] = ppn[ppt==maxdt-2+1]
    nc = nodecost[nodeorder$n[[2,j]]] # get nodecost
    nodeorder$n[[2,j]] = nodeorder$n[[2,j]][order(nc)] #sort by cost
    nodeorder$k[[2,j]] = vector('integer',length(nodeorder$n[[2,j]]))
    nodeorder$k[[2,j]][nodeorder$n[[2,j]] %in% inds] = 1
    #     nodeorder$k[[2,j]] = nodeorder$k[[2,j]][order(nc)]
    
    # For k=0, get ppt==2 nodes
    nodeorder$n[[1,j]] = ppn[ppt==maxdt-1+1]
    nc = nodecost[nodeorder$n[[1,j]]]
    nodeorder$n[[1,j]] = nodeorder$n[[1,j]][order(nc)]
    nodeorder$k[[1,j]] = vector('integer',length(nodeorder$n[[1,j]]))
    #     nodeorder$k[[1,j]][inds] = 1
    #     nodeorder$k[[1,j]] = nodeorder$k[[1,j]][order(nc)]
  }
  return(nodeorder)
}

# --------------------------------------------------------------------------------------------------
#' @title runinf DEPRECATED
#' @description runinf - run inference over all time slices from 1:tend for all specified target 
#' nodes targetnodeord, eliminate all the nodes in that time slice as specified in cptaddch/t (
#' which is nodeorder cpts with repetitions filtered out) using multiply and marginalise. Store
#' the results of inference in the updated joint structure J.<br/>
#' Note that this does FORWARDS INFERENCE only!
#' @param  targetnodeord is a vector containing target nodes $\mathcal{T}$ in the order to be evaled
#' @param tend = index to last t-slice in simulation window 1:tend
#' @param cptaddch[[k=1..maxdt for $\Pi(t-k)$,target node]][[foreach nodeorder node]][vec of 
#'  nodeind corresponding to CPT's to add for multiply step]
#' @param cptaddcht[[k=1..maxdt for $\Pi(t-k)$,target node]][[foreach nodeorder node]][vec of node
#' dt-slice corresponding to CPT's to add for multiply step]
#' @param J[[1..tend,target node]] with $p numeric vec of p, $n int vec of nodes n in $J$ in order
#' in which $J$ is composed, $t int vec of node t-slices 1..tend, $s list of state-labels (as 
#' indices) foreach node in $n$
#' @param Jcpt[[1..tend,target node]] is the original joint populated by $cpt (same struct as J) 
#' @return updated joint structure J (same struct as parameter J)
#' @export
#' @examples none
#' @name runinf
runinf = function(targetnodeord,tend,maxdt,nodeorder,cptaddch,cptaddcht,J,Jcpt) {
  for(i in 1:tend) { # iterate over all t-slices 1...tend
    #   foreach(i=1:tend,.export) %dopar% { # iterate over all t-slices 1...tend
    for(j in 1:length(targetnodeord)) { # iterate over all target nodes
      n = targetnodeord[j] # target node n
      cat('i=',i,'n=',n,'-------------------------------------------------------------------','\n')
      for(t in i:1) { # Backwards inference: foreach t-slice from t to 1
        k = min(i-t+1,maxdt) # index - which t-k setup to use
        if(length(cptaddch[[k,n]])>0) { # If there are nodes to eliminate
          
          for(l in 1:length(cptaddch[[k,n]])) { # foreach node to elim in t-slice
            if(length(cptaddch[[k,n]][[l]])>0) { # If there are nodes to multiply
              
              for(m in 1:length(cptaddch[[k,n]][[l]])) { # foreach node to multiply
                if(t==i && l==1 && m==1) { # For first node in elimination for target node [[i,n]]
                  J[[i,n]] = Jcpt[[t+cptaddcht[[k,n]][[l]][m],cptaddch[[k,n]][[l]][m]]]
                } else {
                  ij = t+cptaddcht[[k,n]][[l]][m]
                  nj = cptaddch[[k,n]][[l]][m]
                  J[[i,n]] = multiply(J[[i,n]]$p,Jcpt[[ij,nj]]$p,J[[i,n]]$n,Jcpt[[ij,nj]]$n,
                                      J[[i,n]]$t,Jcpt[[ij,nj]]$t,J[[i,n]]$s,Jcpt[[ij,nj]]$s)
                }
              } # end foreach node to multiply
            } # end if there are nodes to multiply
            if(length(nodeorder$n[[k,n]])>0)
              J[[i,n]] = marginalise(J[[i,n]]$p,J[[i,n]]$n,J[[i,n]]$t,J[[i,n]]$s,
                                     nodeorder$n[[k,n]][l],t-nodeorder$k[[k,n]][l]) # marginalise
          } # end elim nodes in t-slice
        } else {
          # No nodes to eliminate - do nothing as J is already the marginal Pr
        }
      } # end foreach t-slice in backwards inference
    } # end foreach node
  } # end foreach t-slice
  return(J)
}

# --------------------------------------------------------------------------------------------------
#' @title runlkinf
#' @description run inference over specified target nodes targetnodeord, specified time
#' 1:tend and putting results in joint structure J($p,$n,$t,$s). Like runinf, runlkinf evaluates
#' target nodes going forwards to enable dynamic CPT updating. Unlike runinf, runlkinf divides
#' $\Pi'(t)$ (where $\Pi=\cup_{i=1..t}{\Pi'(i)}$) into subsets/clusters $\pi_j(t)$ where 
#' $\Pi'(i)=\cup_j{\pi_j(t)}$ and $\pi_j\cap\pi_k(t)=\emptyset\forall j\neq k$. Finding the optimal 
#' ordering of nodes AND clusters AND the optimal design of clusters to minimise computation is
#' likely more than NP-hard since optimal ordering is already NP-hard. Also, in order to accommodate
#' forwards-backwards inference, some form of dynamic programming is required. 
#' 
#' Therefore, runlkinf implements the simplest case where there is ONE CLUSTER $pi(t)$ in $\Pi(t)$
#' with link nodes $\lambda(t)$ - i.e. nodes with child(ren) at t+1. Due to bucket/backwards 
#' elimination, $\lambda(t)$ are always eliminated last. This way, instead of having to eliminate 
#' back to t=1, only need to eliminate to t-1. Also, implementation is for maxdt=2 (Markov) system
#' ONLY. Also, linknodes assumed to be sorted in order of elimination (in case there are par-child
#' relationships between link nodes)
#' @param  targetnodeord is a vector containing target nodes $\mathcal{T}$ in the order to be evaled
#' @param tend = index to last t-slice in simulation window 1:tend
#' @param maxdt = (int) max# dt-slices (e.g. Markov network, maxdt = 1; BN, maxdt=0)
#' @param nodeorder - order in which to eliminate nodes; 
#' $n[[k=1..maxdt for $\Pi(t-k)$,target node]] vector of node inds
#' $t[[k=1..maxdt for $\Pi(t-k)$,target node]] vector of node dt-slice t-k
#' @param cptaddch[[k=1..maxdt for $\Pi(t-k)$,target node]][[foreach nodeorder node]][vec of 
#'  nodeind corresponding to CPT's to add for multiply step]
#' @param cptaddcht[[k=1..maxdt for $\Pi(t-k)$,target node]][[foreach nodeorder node]][vec of node
#' dt-slice corresponding to CPT's to add for multiply step]
#' @param J[[1..tend,target node]] with $p numeric vec of p, $n int vec of nodes n in $J$ in order
#' in which $J$ is composed, $t int vec of node t-slices 1..tend, $s list of state-labels (as 
#' indices) foreach node in $n$
#' @param Jcpt[[1..tend,target node]] is the original joint populated by $cpt (same struct as J) 
#' @param Jlink[[1..tend,target node]] is the joint for $\pi(t)$ (same struct as J)
#' @param linknodes = vec of node indices to us as link nodes $\pi(t)$
#' @param lknodeorder - same struct as nodeorder but is for link nodes $\pi(t)$
#' @param lkcptaddch - same struct as cptaddch but is for link nodes $\pi(t)$
#' @param lkcptaddcht - same struct as cptaddcht but for link nodes $\pi(t)$
#' @param updJcpt - function handle for non-homogeneous CPT updating
#' @return list containing updated joint structure J and Jlink ($\pi(t)$)
#' @export
#' @examples none
#' @name runlkinf
runlkinf = function(targetnodeord,tend,maxdt,nodeorder,cptaddch,cptaddcht,J,Jcpt,
                    Jlink,linknodes,lknodeorder,lkcptaddch,lkcptaddcht,updJcpt) {
  # Perform backwards elimination inference - first do Jlink
  Jlink = vector('list',tend) # Initialise Jlink struct foreach t=1..tend
  Jlink = lapply(Jlink,function(X) list(p=numeric(0),n=integer(0),t=integer(0),s=integer(0)))
  et = proc.time()[3]
  for(i in 1:tend) { # iterate over all t-slices 1...tend
    if(i != 1 && !is.null(updJcpt)) 
      Jcpt = updJcpt(i,J,Jcpt,cptaddch,linknodes,Jlink) # Dynamically update Jcpt using custom fcn
    cat('i=',i,'-----------------------------------------\n')
    
    # Eliminate nodes for linknodes and store to dedicated Jlink
    for(l in 1:length(lkcptaddch)) { # foreach node to elim in t-slice
      if(length(lkcptaddch[[l]])>0) { 
        for(m in 1:length(lkcptaddch[[l]])) { # foreach node to multiply
          if(l==1 && m==1) { # For first node in elimination
            Jlink[[i]] = Jcpt[[i+lkcptaddcht[[l]][m],lkcptaddch[[l]][m]]]
          } else { 
            ij = i+lkcptaddcht[[l]][m]
            nj = lkcptaddch[[l]][m]
            Jlink[[i]] = multiply(Jlink[[i]]$p,Jcpt[[ij,nj]]$p,Jlink[[i]]$n,Jcpt[[ij,nj]]$n,
                                  Jlink[[i]]$t,Jcpt[[ij,nj]]$t,Jlink[[i]]$s,Jcpt[[ij,nj]]$s)
          }
        } # end foreach node to multiply
      } # end if there are nodes to multiply
      # Assumes there is always a node to eliminate
      Jlink[[i]] = marginalise(Jlink[[i]]$p,Jlink[[i]]$n,Jlink[[i]]$t,Jlink[[i]]$s,
                               lknodeorder$n[[1]][l],i-lknodeorder$k[[1]][l]) # marginalise
    } # end foreach nodeorder node
    
    # Eliminate \lambda link nodes - these can be done at the end due to backwards elimination
    if(i != 1) {
      # Multiply joint by link nodes J from t-1
      Jlink[[i]] = multiply(Jlink[[i-1]]$p,Jlink[[i]]$p,Jlink[[i-1]]$n,Jlink[[i]]$n,
                            Jlink[[i-1]]$t,Jlink[[i]]$t,Jlink[[i-1]]$s,Jlink[[i]]$s)
      for(l in 1:length(linknodes)) { # Marginalise by link nodes from t-1
        Jlink[[i]] = marginalise(Jlink[[i]]$p,Jlink[[i]]$n,Jlink[[i]]$t,Jlink[[i]]$s,
                                 linknodes[l],i-1) # marginalise
      }
    } # end eliminating \lambda link nodes from t-1
    
    # Perform backwards elimination inference - now do target nodes at each t-slice
    for(j in 1:length(targetnodeord)) { # iterate over all target nodes
      n = targetnodeord[j]
      k = 1 # index - which t-k setup to use - always use 1 for maxdt=2 (Markov system)
      if(length(cptaddch[[k,n]])>0) { # If there are nodes to eliminate
        for(l in 1:length(cptaddch[[k,n]])) { # foreach node to elim in t-slice t
          if(length(cptaddch[[k,n]][[l]])>0) { # If there are nodes to multiply
            for(m in 1:length(cptaddch[[k,n]][[l]])) { # foreach node to multiply
              if(l==1 && m==1) { # For first node in elimination for target node, add CPT
                J[[i,n]] = Jcpt[[i+cptaddcht[[k,n]][[l]][m],cptaddch[[k,n]][[l]][m]]]
              } else {
                ij = i + cptaddcht[[k,n]][[l]][m]
                nj = cptaddch[[k,n]][[l]][m]
                J[[i,n]] = multiply(J[[i,n]]$p,Jcpt[[ij,nj]]$p,J[[i,n]]$n,Jcpt[[ij,nj]]$n,
                                    J[[i,n]]$t,Jcpt[[ij,nj]]$t,J[[i,n]]$s,Jcpt[[ij,nj]]$s)
              }
            } # end foreach node to multiply
          } # end if there are nodes to multiply
          if(length(nodeorder$n[[k,n]])>0)
            J[[i,n]] = marginalise(J[[i,n]]$p,J[[i,n]]$n,J[[i,n]]$t,J[[i,n]]$s,
                                   nodeorder$n[[k,n]][l],i-nodeorder$k[[k,n]][l]) # marginalise
        } # end foreach node to elim
        # Eliminate \lambda link nodes - these can be done at the end due to backwards elimination
        if(i != 1) {
          # Multiply joint by link nodes J from t-1
          J[[i,n]] = multiply(Jlink[[i-1]]$p,J[[i,n]]$p,Jlink[[i-1]]$n,J[[i,n]]$n,
                              Jlink[[i-1]]$t,J[[i,n]]$t,Jlink[[i-1]]$s,J[[i,n]]$s)
          for(l in 1:length(linknodes)) { # Marginalise by link nodes from t-1
            J[[i,n]] = marginalise(J[[i,n]]$p,J[[i,n]]$n,J[[i,n]]$t,J[[i,n]]$s,
                                   linknodes[l],i-1) # marginalise
          }
        } # end eliminating \lambda link nodes from t-1
      } else {
        J[[i,n]] = Jcpt[[i,n]] # No nodes to eliminate - CPT is already the marginal Pr
      } # end if there are nodes to elim
    } # end foreach target node
  } # end foreach t-slice
  cat('Inference for target nodes took  ',proc.time()[3]-et,'seconds\n')
  return(list(J=J,Jlink=Jlink,Jcpt=Jcpt))
}

#' @title setupPi
#' @description Create the custom Pinodeorder, Picptaddch, Picptaddcht structures to be used in 
#' oneslicemargs. 
#' @param DT - dt-slice index to use for childnodeind
#' @param nodecost[1..maxdt,node] containing cost value $OC=D+C$
#' @param childnodeind ([[1..maxdt,node]][child node ind to nnames])
#' @param fbmode - 0 for forwards only, 1 for forwards and backwards
#' @return list containing: $Pinodeorder = is equivalent to lknodeorder wrapped by layer d, 
#' Pinodeorder[[d]]= defines the order of elimination of nodes in this layer given results from prev
#' layers; $n = node id, $k = dt-slice. $Picptaddch[[layer d]]  contains list of vec of nodeinds
#' corresponding to CPT's to add for multiply foreach elimination step in Pibknodeorder[[d]].
#' $Picptaddcht[[layer d]]- ditto but for node dt-slice
#' @export
#' @examples none
#' @name setupPi
setupPi = function(DT,nodecost,childnodeind,fbmode) {
  ordnodecost = order(nodecost[2,])
  intnodecost = floor(nodecost[2,ordnodecost])
  
  # ASSUMES: $\Pi(t)$ so finding $P(X(t)|\pmb{E})$ all nodes are in same time slice so all $k/t ==0
  # Propagate wave backwards foreach $\pi_d(t)$ from d=0 to d=D-1 BUT since R indexing starts at 1,
  # translate [0,D-1] to [1,D]. 
  Pinodeorder = NULL
  Picptaddch = NULL
  Picptaddcht = NULL
  # For forwards and backwards inference
  for(i in 1:length(unique(floor(nodecost[2,]))) ) { # foreach d-layer i (target layer for inference)
    Pinodeorder[[i]] = list()
    if(fbmode==0) {
      Pinodeorder[[i]]$n = ordnodecost[intnodecost>i-1] # Elim all nodes in all layers except layer d
    } else if(fbmode==1) {
      Pinodeorder[[i]]$n = ordnodecost[intnodecost!=i-1] # Elim all nodes in all layers except layer d
    }
    Pinodeorder[[i]]$k = 0*Pinodeorder[[i]]$n
    # Because of backwards elimination, we first assume that all nodes at future slices have already
    # been eliminated. Also, all child nodes are eliminated first (by layer starting with layer 1
    Picptaddch[[i]] = sapply(Pinodeorder[[i]]$n,function(X)list(X))
    nd = ordnodecost[intnodecost==i-1] # nodes to add in layer d
    # For layer d+1, need to also multiply by nodes in layer d (target layer for inference)
    if(i<length(unique(floor(nodecost[2,])))) {
      ind1 = which(Pinodeorder[[i]]$n==ordnodecost[intnodecost==i][1]) # Get indices to layer d+1
      ind2=which(Pinodeorder[[i]]$n==ordnodecost[intnodecost==i][length(ordnodecost[intnodecost==i])])
      for(j in ind1:ind2) {
        prev = unlist(Picptaddch[[i]][1:j])
        if(j<ind2) {
          Picptaddch[[i]][[j]]=c(Picptaddch[[i]][[j]],
                                 setdiff(intersect(childnodeind[[DT,Pinodeorder[[i]]$n[j]]],nd),prev))
        } else {
          # In the last element, add all remaining nodes from layer d to cptaddch for d+1
          Picptaddch[[i]][[j]]=c(Picptaddch[[i]][[j]],setdiff(nd,prev))
        }
      }
    } else { # layer d=D
      ind = length(Picptaddch[[i]]) # ind at minimum is 1
      if(ind>0) {
        Picptaddch[[i]][[ind]] = c(Picptaddch[[i]][[ind]],nd)
      } else {
        Picptaddch[[i]] = list(nd)
      }
    }
    Picptaddcht[[i]] = sapply(Picptaddch[[i]],function(X)list(X*0))
  }
  return(list(Pinodeorder=Pinodeorder,Picptaddch=Picptaddch,Picptaddcht=Picptaddcht))
}


