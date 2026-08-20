amerge <- function(d1, d2, firm=NULL, approx=NULL) {
  rt = Sys.time()
  
  # Take care of conflicting column names
  n2 = data.frame(oldname = names(d2), newname = names(d2))
  n2$newname = as.character(n2$newname)
  n2$newname[(n2$oldname %in% names(d1)) & !(n2$oldname %in% firm)] =
    paste(n2$newname[(n2$oldname %in% names(d1)) & !(n2$oldname %in% firm)], "2", sep=".")
  
  # Add unique row IDs
  if (length(c(firm, approx))>1) {
    d1$ID1 = factor(apply(d1[,c(approx,firm)], 1, paste, collapse=" "))
    d2$ID2 = factor(apply(d2[,c(approx,firm)], 1, paste, collapse=" "))
  } else {
    d1$ID1 = factor(d1[,c(approx,firm)])
    d2$ID2 = factor(d2[,c(approx,firm)])
  }
  
  # Perform initial merge on the 'firm' parameters, if any
  # Otherwise match all to all
  if (length(firm)>0) {
    t1 = merge(d1, d2, by=firm, all=T, suff=c("",".2"))
  } else {
    names(d2)= c(n2$newname,"ID2")
    t1 = data.frame()
    for (i1 in 1:nrow(d1)) {
      trow = d1[i1,]
      t1 = rbind(t1, cbind(trow, d2))
    }
  }
  
  # Match by the most approximate record
  if (length(approx)==1) {
    # Calculate the differential for approximate merging
    t1$DIFF = abs(t1[,approx] - t1[,n2$newname[n2$oldname==approx]])
    # Sort data by ascending DIFF, so that best matching records are used first
    t1 = t1[order(t1$DIFF, t1$ID1, t1$ID2),]
    t2 = data.frame()
    d2$used = 0
    # For each record of d1, find match from d2
    for (i1 in na.omit(unique(t1$ID1))) {
      tx = t1[!is.na(t1$DIFF) & t1$ID1==i1,]
      # If there are non-missing records, get the one with minimum DIFF (top one)
      if (nrow(tx)>0) {
        tx = tx[1,]
        # If matching record found, remove it from the pool, so it's not used again
        t1[!is.na(t1$ID2) & t1$ID2==tx$ID2, c(n2$newname[!(n2$newname %in% firm)], "DIFF")] = NA
        # And mark it as used
        d2$used[d2$ID2==tx$ID2] = 1
      } else {
        # If there are no non-missing records, just get the first one from the top
        tx = t1[!is.na(t1$ID1) & t1$ID1==i1,][1,]
      }
      t2 = rbind(t2,tx)
    }
  } else {
    t2 = t1
  }
  # Make the records the same order as d1
  t2 = t2[match(d1$ID1, t2$ID1),]
  # Add unmatched records from d2 to the end of output
  if (any(d2$used==0)) {
    tx = t1[t1$ID2 %in% d2$ID2[d2$used==0], ]
    tx = tx[!duplicated(tx$ID2),]
    tx[, names(d1)[!(names(d1) %in% c(firm))]] = NA
    t2 = rbind(t2,tx)
    t2[is.na(t2[,approx]), approx] = t2[is.na(t2[,approx]), n2$newname[n2$oldname==approx]]
  }
  t2$DIFF = t2$ID1 = t2$ID2 = NULL
  cat("* Run time: ", round(difftime(Sys.time(),rt, "secs"),1), " seconds.\n", sep="")
  return(t2)
}