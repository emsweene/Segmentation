

quick_auc = function(pred, 
  fpr.stop = c(0.005, 0.01, 1)
  ) {

  xvals = pred@fp[[1]] / pred@n.neg[[1]]
  yvals = pred@tp[[1]] / pred@n.pos[[1]]
  finite.bool <- is.finite(xvals) & is.finite(yvals)
  xvals = xvals[finite.bool]
  yvals = yvals[finite.bool]
  cum_add = function(y) {
    n = length(y)
    y[2:n] + y[1:(n-1)]
  }  
  aucs = lapply(fpr.stop, function(r) {
    ind = max(which(xvals <= r))
    if (ind == length(xvals)) {
      tpr.stop = yvals[ind]
    } else {
    tpr.stop <- approxfun(
      xvals[ind:(ind + 1)], 
      yvals[ind:(ind + 1)])(r)
    }
    xx <- c(xvals[1:ind], r)
    yy <- c(yvals[1:ind], tpr.stop)
    auc = sum(diff(xx) * cum_add(yy) * 0.5)
    auc = auc / r
    return(auc)
  })
  names(aucs) = fpr.stop
  print(aucs)
  return(aucs)
}


pred_binary = function(predictions, 
  labels, check_lab = TRUE) {
  npred = NROW(predictions)
  nlab = NROW(labels)
  if (npred != nlab) {
    stop("Numbers don't add up - data mismatch")
  }
  if (length(labels) != nlab ||
    !is.atomic(labels)) {
    stop("it doesn't seem that the labels are a vector")
  }
  if (check_lab) {
    ulab = unique(labels)
    ulab = as.numeric(ulab)
    if (!all(ulab %in% c(0, 1))) {
      stop(paste0("Either missing data or non 0/1 data",
        " exist - must fix on your own"))
    }
  }
  labels = as.logical(labels)
  n.pos = sum(labels)
  n.neg = sum(1 - labels)
  ans = fast_roc(xdata$oasis_p, labels = yy)

  fast_roc = function(
    predictions, 
    labels,
    check_lab = FALSE) {

    # fastest sort - decreasing
    df = data.table::data.table(
      x = -predictions,
      y = labels,
      key = "x")
    df = dplyr::as.tbl(df)
    df$x = -df$x

    tp = cumsum(df$y)
    fp = 1:nrow(df) - tp
    notdups <- !rev(duplicated(rev(df$x)))
    tp <- c(0, tp[notdups])
    fp <- c(0, fp[notdups])
    cutoffs <- c(Inf, df$x[notdups])
    return(list(
      cutoffs = cutoffs, fp = fp, tp = tp,
      n.pos.pred = tp + fp))
  }

  npred = ncol(predictions)
  seq_n = seq(npred)
  n.pos = lapply(seq_n, function(r) n.pos)
  n.neg = lapply(seq_n, function(r) n.neg)

  fp = tp = vector(
    mode = "list", 
    length = npred)
  n.pos.pred = cutoffs = tn = fp
  n.neg.pred = fn = fp

  for (i in 1:length(predictions)) {
      ans <- fast_roc(predictions[[i]], 
          labels[[i]])
      cutoffs <- c(cutoffs, list(ans$cutoffs))
      fp <- c(fp, list(ans$fp))
      tp <- c(tp, list(ans$tp))
      fn <- c(fn, list(n.pos[[i]] - tp[[i]]))
      tn <- c(tn, list(n.neg[[i]] - fp[[i]]))
      n.pos.pred <- c(n.pos.pred, list(tp[[i]] + fp[[i]]))
      n.neg.pred <- c(n.neg.pred, list(tn[[i]] + fn[[i]]))
  }

  new("prediction", 
    predictions = predictions, 
    labels = labels, 
    cutoffs = cutoffs, 
    fp = fp, tp = tp, 
    fn = fn, tn = tn, 
    n.pos = n.pos, n.neg = n.neg, n.pos.pred = n.pos.pred, 
    n.neg.pred = n.neg.pred)
}