###### -- Pseudogene Congruence between assemblers ----------------------------

###### -- libraries -----------------------------------------------------------

suppressMessages(library(SynExtend))
suppressMessages(library(VennDiagram))

###### -- ad hoc functions ----------------------------------------------------

# specific adjustments to VennDiagram's major function for adding a 4 option venn
# diagram to another plot...
venn.diagram.adhoc <- function(
    x,
    filename,
    disable.logging = FALSE,
    height = 3000,
    width = 3000,
    resolution = 500,
    imagetype = 'tiff',
    units = 'px',
    compression = 'lzw',
    na = 'stop',
    main = NULL,
    sub = NULL,
    main.pos = c(0.5, 1.05),
    main.fontface = 'plain',
    main.fontfamily = 'serif',
    main.col = 'black',
    main.cex = 1,
    main.just = c(0.5, 1),
    sub.pos = c(0.5, 1.05),
    sub.fontface = 'plain',
    sub.fontfamily = 'serif',
    sub.col = 'black',
    sub.cex = 1,
    sub.just = c(0.5, 1),
    category.names = names(x),
    force.unique = TRUE,
    print.mode = 'raw',
    sigdigs = 3,
    direct.area = FALSE,
    area.vector = 0,
    hyper.test = FALSE,
    total.population = NULL,
    lower.tail = TRUE,
    N_Sets = NULL,
    ...
) {
  
  #Create a string to capture the date and the time of day
  time.string = gsub(':', '-', gsub(' ', '_', as.character(Sys.time())))
  
  #Initialize the logger to output to file
  if (disable.logging) {
    flog.appender(appender.console(), name = 'VennDiagramLogger');
  } else {
    flog.appender(appender.file(
      paste0(if (!is.null(filename)) filename else 'VennDiagram', '.', time.string, '.log')),
      name = 'VennDiagramLogger'
    );
  }
  
  #Log the parameters the function was called with
  out.list = as.list(sys.call())
  out.list[[1]] <- NULL
  out.string = capture.output(out.list)
  flog.info(out.string,name='VennDiagramLogger')
  
  #If the input area.vector correspond directly to a1,a2, etc, then call the function and pass it through directly by a flag
  if(direct.area){
    if(1 == length(area.vector)){
      list.names <- category.names;
      if (is.null(list.names)) { list.names <- ''; }
      grob.list <- VennDiagram::draw.single.venn(
        area = area.vector[1],
        category = list.names,
        ind = FALSE,
        ...
      );
    }
    if(3 == length(area.vector)){
      grob.list <- VennDiagram::draw.pairwise.venn(
        area1 = area.vector[1],
        area2 = area.vector[2],
        cross.area = area.vector[3],
        category = category.names,
        ind = FALSE,
        print.mode=print.mode,
        sigdigs=sigdigs,
        ...
      );
    }
    if(7 == length(area.vector)){
      grob.list <- VennDiagram::draw.triple.venn(
        area1 = 0,
        area2 = 0,
        area3 = 0,
        n12 = 0,
        n23 = 0,
        n13 = 0,
        n123 = 0,
        category = category.names,
        ind = FALSE,
        list.order = 1:3,
        print.mode=print.mode,
        sigdigs=sigdigs,
        area.vector=area.vector,
        direct.area=TRUE,
        ...
      );
    }
    if(15 == length(area.vector)){
      grob.list <- VennDiagram::draw.quad.venn(
        area1 = 0,
        area2 = 0,
        area3 = 0,
        area4 = 0,
        n12 = 0,
        n13 = 0,
        n14 = 0,
        n23 = 0,
        n24 = 0,
        n34 = 0,
        n123 = 0,
        n124 = 0,
        n134 = 0,
        n234 = 0,
        n1234 = 0,
        category = category.names,
        ind = FALSE,
        print.mode=print.mode,
        sigdigs=sigdigs,
        area.vector=area.vector,
        direct.area=TRUE,
        ...
      );
    }
    if(31 == length(area.vector)){
      grob.list <- VennDiagram::draw.quintuple.venn(
        area1 = 0,
        area2 = 0,
        area3 = 0,
        area4 = 0,
        area5 = 0,
        n12 = 0,
        n13 = 0,
        n14 = 0,
        n15 = 0,
        n23 = 0,
        n24 = 0,
        n25 = 0,
        n34 = 0,
        n35 = 0,
        n45 = 0,
        n123 = 0,
        n124 = 0,
        n125 = 0,
        n134 = 0,
        n135 = 0,
        n145 = 0,
        n234 = 0,
        n235 = 0,
        n245 = 0,
        n345 = 0,
        n1234 = 0,
        n1235 = 0,
        n1245 = 0,
        n1345 = 0,
        n2345 = 0,
        n12345 = 0,
        category = category.names,
        ind = FALSE,
        print.mode=print.mode,
        sigdigs=sigdigs,
        area.vector=area.vector,
        direct.area=TRUE,
        ...
      );
    }
  }
  
  #Use default processing behaviour of having individual elements in the list x
  else{
    
    if (force.unique) {
      for (i in 1:length(x)) {
        x[[i]] <- unique(x[[i]])
      }
    }
    
    # check for the presence of NAs in the input list
    if ('none' == na) {
      x <- x;
    }
    else if ('stop' == na) {
      for (i in 1:length(x)) {
        # stop if there are any NAs in this vector
        if (any(is.na(x[[i]]))) {
          flog.error('NAs in dataset', call. = FALSE,name='VennDiagramLogger')
          stop('NAs in dataset', call. = FALSE);
        }
      }
    }
    else if ('remove' == na) {
      for (i in 1:length(x)) { x[[i]] <- x[[i]][!is.na(x[[i]])]; }
    }
    else {
      flog.error('Invalid na option: valid options are "none", "stop", and "remove"',name='VennDiagramLogger')
      stop('Invalid na option: valid options are "none", "stop", and "remove"');
    }
    
    # check the length of the given list
    if (0 == length(x) | length(x) > 5) {
      flog.error('Incorrect number of elements.', call. = FALSE,name='VennDiagramLogger')
      stop('Incorrect number of elements.', call. = FALSE);
    }
    
    # draw a single-set Venn diagram
    if (1 == length(x)) {
      list.names <- category.names;
      if (is.null(list.names)) { list.names <- ''; }
      grob.list <- VennDiagram::draw.single.venn(
        area = length(x[[1]]),
        category = list.names,
        ind = FALSE,
        ...
      );
    }
    
    # draw a pairwise Venn diagram
    else if (2 == length(x)) {
      grob.list <- VennDiagram::draw.pairwise.venn(
        area1 = length(x[[1]]),
        area2 = length(x[[2]]),
        cross.area = length(intersect(x[[1]],x[[2]])),
        category = category.names,
        ind = FALSE,
        print.mode=print.mode,
        sigdigs=sigdigs,
        ...
      );
    }
    
    # draw a three-set Venn diagram
    else if (3 == length(x)) {
      A <- x[[1]];
      B <- x[[2]];
      C <- x[[3]];
      
      list.names <- category.names;
      
      nab <- intersect(A, B);
      nbc <- intersect(B, C);
      nac <- intersect(A, C);
      
      nabc <- intersect(nab, C);
      
      grob.list <- VennDiagram::draw.triple.venn(
        area1 = length(A),
        area2 = length(B),
        area3 = length(C),
        n12 = length(nab),
        n23 = length(nbc),
        n13 = length(nac),
        n123 = length(nabc),
        category = list.names,
        ind = FALSE,
        list.order = 1:3,
        print.mode=print.mode,
        sigdigs=sigdigs,
        ...
      );
    }
    
    # draw a four-set Venn diagram
    else if (4 == length(x)) {
      A <- x[[1]];
      B <- x[[2]];
      C <- x[[3]];
      D <- x[[4]];
      
      list.names <- category.names;
      
      n12 <- intersect(A, B);
      n13 <- intersect(A, C);
      n14 <- intersect(A, D);
      n23 <- intersect(B, C);
      n24 <- intersect(B, D);
      n34 <- intersect(C, D);
      
      n123 <- intersect(n12, C);
      n124 <- intersect(n12, D);
      n134 <- intersect(n13, D);
      n234 <- intersect(n23, D);
      
      n1234 <- intersect(n123, D);
      
      grob.list <- draw.quad.venn(
        area1 = length(A),
        area2 = length(B),
        area3 = length(C),
        area4 = length(D),
        n12 = length(n12),
        n13 = length(n13),
        n14 = length(n14),
        n23 = length(n23),
        n24 = length(n24),
        n34 = length(n34),
        n123 = length(n123),
        n124 = length(n124),
        n134 = length(n134),
        n234 = length(n234),
        n1234 = length(n1234),
        category = list.names,
        ind = FALSE,
        print.mode=print.mode,
        sigdigs=sigdigs,
        N_Sets = N_Sets,
        ...
      );
    }
    
    # draw a five-set Venn diagram
    else if (5 == length(x)) {
      A <- x[[1]];
      B <- x[[2]];
      C <- x[[3]];
      D <- x[[4]];
      E <- x[[5]];
      
      list.names <- category.names;
      
      n12 <- intersect(A, B);
      n13 <- intersect(A, C);
      n14 <- intersect(A, D);
      n15 <- intersect(A, E);
      n23 <- intersect(B, C);
      n24 <- intersect(B, D);
      n25 <- intersect(B, E);
      n34 <- intersect(C, D);
      n35 <- intersect(C, E);
      n45 <- intersect(D, E);
      
      n123 <- intersect(n12, C);
      n124 <- intersect(n12, D);
      n125 <- intersect(n12, E);
      n134 <- intersect(n13, D);
      n135 <- intersect(n13, E);
      n145 <- intersect(n14, E);
      n234 <- intersect(n23, D);
      n235 <- intersect(n23, E);
      n245 <- intersect(n24, E);
      n345 <- intersect(n34, E);
      
      n1234 <- intersect(n123, D);
      n1235 <- intersect(n123, E);
      n1245 <- intersect(n124, E);
      n1345 <- intersect(n134, E);
      n2345 <- intersect(n234, E);
      
      n12345 <- intersect(n1234, E);
      
      grob.list <- VennDiagram::draw.quintuple.venn(
        area1 = length(A),
        area2 = length(B),
        area3 = length(C),
        area4 = length(D),
        area5 = length(E),
        n12 = length(n12),
        n13 = length(n13),
        n14 = length(n14),
        n15 = length(n15),
        n23 = length(n23),
        n24 = length(n24),
        n25 = length(n25),
        n34 = length(n34),
        n35 = length(n35),
        n45 = length(n45),
        n123 = length(n123),
        n124 = length(n124),
        n125 = length(n125),
        n134 = length(n134),
        n135 = length(n135),
        n145 = length(n145),
        n234 = length(n234),
        n235 = length(n235),
        n245 = length(n245),
        n345 = length(n345),
        n1234 = length(n1234),
        n1235 = length(n1235),
        n1245 = length(n1245),
        n1345 = length(n1345),
        n2345 = length(n2345),
        n12345 = length(n12345),
        category = list.names,
        ind = FALSE,
        print.mode=print.mode,
        sigdigs=sigdigs,
        ...
      );
    }
    
    # this should never happen because of the previous check
    else {
      flog.error('Invalid size of input object',name='VennDiagramLogger')
      stop('Invalid size of input object');
    }
  }
  
  # if there are two sets in the VennDiagram and the hypergeometric test is requested then perform the test and add the pvalue to the subtitle
  #p value always shown with 2 sig digs. Add another parameter for this later if you want to control the sig digs
  
  if (length(x) == 2 & !is.null(total.population) & hyper.test){
    val.p = calculate.overlap.and.pvalue(x[[1]],x[[2]],total.population, lower.tail = lower.tail);
    if(is.null(sub)){
      sub = paste0('p = ',signif(val.p[3],digits=2))
    }else{
      sub = paste0(sub,', p = ',signif(val.p[3],digits=2))
    }
  }
  
  # if requested, add a sub-title
  if (!is.null(sub)) {
    grob.list <- add.title(
      gList = grob.list,
      x = sub,
      pos = sub.pos,
      fontface = sub.fontface,
      fontfamily = sub.fontfamily,
      col = sub.col,
      cex = sub.cex
    );
  }
  
  # if requested, add a main-title
  if (!is.null(main)) {
    grob.list <- add.title(
      gList = grob.list,
      x = main,
      pos = main.pos,
      fontface = main.fontface,
      fontfamily = main.fontfamily,
      col = main.col,
      cex = main.cex
    );
  }
  
  # just return the plot into the graphics device...
  grid.draw(grob.list)
  
  # return(c(length(A),
  #          length(B),
  #          length(C),
  #          length(D),
  #          length(n12),
  #          length(n13),
  #          length(n14),
  #          length(n23),
  #          length(n24),
  #          length(n34),
  #          length(n123),
  #          length(n124),
  #          length(n134),
  #          length(n234),
  #          length(n1234)))
  
  # if a filename is given, write a desired image type, TIFF default
  # if (!is.null(filename)) {
  #   
  #   # set the graphics driver
  #   current.type <- getOption('bitmapType');
  #   if (length(grep('Darwin', Sys.info()['sysname']))) {
  #     options(bitmapType = 'quartz');
  #   }
  #   else {
  #     options(bitmapType = 'cairo');
  #   }
  #   
  #   # TIFF image type specified
  #   if('tiff' == imagetype) {
  #     tiff(
  #       filename = filename,
  #       height = height,
  #       width = width,
  #       units = units,
  #       res = resolution,
  #       compression = compression
  #     );
  #   }
  #   
  #   # PNG image type specified
  #   else if('png' == imagetype) {
  #     png(
  #       filename = filename,
  #       height = height,
  #       width = width,
  #       units = units,
  #       res = resolution
  #     );
  #   }
  #   
  #   # SVG image type specified
  #   else if('svg' == imagetype) {
  #     svg(
  #       filename = filename,
  #       height = height,
  #       width = width
  #     );
  #   }
  #   
  #   # Invalid imagetype specified
  #   else {
  #     flog.error('You have misspelled your "imagetype", please try again',name='VennDiagramLogger')
  #     stop('You have misspelled your "imagetype", please try again');
  #   }
  #   
  #   grid.draw(grob.list);
  #   dev.off();
  #   options(bitmapType = current.type);
  #   
  #   # return a success code
  #   return(1);
  # }
  # 
  # # if file creation was not requested return the plotting object
  # return(grob.list);
}

### FUNCTION TO DRAW VENN DIAGRAM WITH FOUR SETS #################################################
draw.quad.venn <- function(
    area1,
    area2,
    area3,
    area4,
    n12,
    n13,
    n14,
    n23,
    n24,
    n34,
    n123,
    n124,
    n134,
    n234,
    n1234,
    category = rep('', 4),
    lwd = rep(2, 4),
    lty = rep('solid', 4),
    col = rep('black', 4),
    fill = NULL,
    alpha = rep(0.5, 4),
    label.col = rep('black', 15),
    cex = rep(1, 15),
    fontface = rep('plain', 15),
    fontfamily = rep('serif', 15),
    cat.pos = c(-15, 15, 0, 0),
    cat.dist = c(0.22, 0.22, 0.11, 0.11),
    cat.col = rep('black', 4),
    cat.cex = rep(1, 4),
    cat.fontface = rep('plain', 4),
    cat.fontfamily = rep('serif', 4),
    cat.just = rep(list(c(0.5, 0.5)), 4),
    rotation.degree = 0,
    rotation.centre = c(0.5, 0.5),
    ind = TRUE,
    cex.prop=NULL,
    print.mode = 'raw',
    sigdigs=3,
    direct.area = FALSE,
    area.vector = 0,
    N_Sets = NULL,
    ...
) {
  
  #area1 > area2 > area3 > area4
  # check parameter lengths
  if (length(category) == 1) { cat <- rep(category, 4); }
  else if (length(category) != 4) { flog.error('Unexpected parameter length for "category"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "category"'); }
  
  if (length(lwd) == 1) { lwd <- rep(lwd, 4); }
  else if (length(lwd) != 4) { flog.error('Unexpected parameter length for "lwd"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "lwd"'); }
  
  if (length(lty) == 1) { lty <- rep(lty, 4); }
  else if (length(lty) != 4) { flog.error('Unexpected parameter length for "lty"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "lty"'); }
  
  if (length(col) == 1) { col <- rep(col, 4); }
  else if (length(col) != 4) { flog.error('Unexpected parameter length for "col"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "col"'); }
  
  if (length(label.col) == 1) { label.col <- rep(label.col, 15); }
  else if (length(label.col) != 15) { flog.error('Unexpected parameter length for "label.col"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "label.col"'); }
  
  if (length(cex) == 1) { cex <- rep(cex, 15); }
  else if (length(cex) != 15) { flog.error('Unexpected parameter length for "cex"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cex"'); }
  
  if (length(fontface) == 1) { fontface <- rep(fontface, 15); }
  else if (length(fontface) != 15) { flog.error('Unexpected parameter length for "fontface"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "fontface"'); }
  
  if (length(fontfamily) == 1) { fontfamily <- rep(fontfamily, 15); }
  else if (length(fontfamily) != 15) { flog.error('Unexpected parameter length for "fontfamily"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "fontfamily"'); }
  
  if (length(fill) == 1) { fill <- rep(fill, 4); }
  else if (length(fill) != 4 & length(fill) != 0) { flog.error('Unexpected parameter length for "fill"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "fill"'); }
  
  if (length(alpha) == 1) { alpha <- rep(alpha, 4); }
  else if (length(alpha) != 4 & length(alpha) != 0) { flog.error('Unexpected parameter length for "alpha"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "alpha"'); }
  
  if (length(cat.pos) == 1) { cat.pos <- rep(cat.pos, 4); }
  else if (length(cat.pos) != 4) { flog.error('Unexpected parameter length for "cat.pos"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.pos"'); }
  
  if (length(cat.dist) == 1) { cat.dist <- rep(cat.dist, 4); }
  else if (length(cat.dist) != 4) { flog.error('Unexpected parameter length for "cat.dist"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.dist"'); }
  
  if (length(cat.col) == 1) { cat.col <- rep(cat.col, 4); }
  else if (length(cat.col) != 4) { flog.error('Unexpected parameter length for "cat.col"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.col"'); }
  
  if (length(cat.cex) == 1) { cat.cex <- rep(cat.cex, 4); }
  else if (length(cat.cex) != 4) { flog.error('Unexpected parameter length for "cat.cex"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.cex"'); }
  
  if (length(cat.fontface) == 1) { cat.fontface <- rep(cat.fontface, 4); }
  else if (length(cat.fontface) != 4) { flog.error('Unexpected parameter length for "cat.fontface"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.fontface"'); }
  
  if (length(cat.fontfamily) == 1) { cat.fontfamily <- rep(cat.fontfamily, 4); }
  else if (length(cat.fontfamily) != 4) { flog.error('Unexpected parameter length for "cat.fontfamily"',name='VennDiagramLogger')
    stop('Unexpected parameter length for "cat.fontfamily"'); }
  
  if (!(is.list(cat.just) && length(cat.just) == 4 && length(cat.just[[1]]) == 2 && length(cat.just[[2]]) == 2 && length(cat.just[[3]]) == 2 && length(cat.just[[4]]) == 2)) { flog.error('Unexpected parameter format for "cat.just"',name='VennDiagramLogger')
    stop('Unexpected parameter format for "cat.just"'); }
  cat.pos <- cat.pos + rotation.degree;
  
  if(direct.area){
    areas <- area.vector;
    #create the variables and assign their values from the area vector
    for(i in 1:15)
    {
      assign(paste('a',i,sep=''),area.vector[i]);
    }
  }
  else {
    # generate partial areas from given arguments
    a6  <- n1234;
    a12 <- n123 - a6;
    a11 <- n124 - a6;
    a5  <- n134 - a6;
    a7  <- n234 - a6;
    a15 <- n12 - a6 - a11 - a12;
    a4  <- n13 - a6 - a5 - a12;
    a10 <- n14 - a6 - a5 - a11;
    a13 <- n23 - a6 - a7 - a12;
    a8  <- n24 - a6 - a7 - a11;
    a2  <- n34 - a6 - a5 - a7;
    a9  <- area1 - a4 - a5 - a6 - a10 - a11 - a12 - a15;
    a14 <- area2 - a6 - a7 - a8 - a11 - a12 - a13 - a15;
    a1  <- area3 - a2 - a4 - a5 - a6 - a7 - a12 - a13;
    a3  <- area4 - a2 - a5 - a6 - a7 - a8 - a10 - a11;
    
    # check plausibility and 0 partial areas
    # areas <- c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15);
    areas <- signif(c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15) / N_Sets,
                    digits = 4L)
  }
  areas.error <- c(
    'a1  <- area3 - a2 - a4 - a5 - a6 - a7 - a12 - a13',
    'a2  <- n34 - a6 - a5 - a7',
    'a3  <- area4 - a2 - a5 - a6 - a7 - a8 - a10 - a11',
    'a4  <- n13 - a6 - a5 - a12',
    'a5  <- n134 - a6',
    'a6  <- n1234',
    'a7  <- n234 - a6',
    'a8  <- n24 - a6 - a7 - a11',
    'a9  <- area1 - a4 - a5 - a6 - a10 - a11 - a12 - a15',
    'a10 <- n14 - a6 - a5 - a11',
    'a11 <- n124 - a6',
    'a12 <- n123 - a6',
    'a15 <- n12 - a6 - a11 - a12',
    'a13 <- n23 - a6 - a7 - a12',
    'a14 <- area2 - a6 - a7 - a8 - a11 - a12 - a13 - a15'
  );
  for (i in 1:length(areas)) {
    if (areas[i] < 0) {
      flog.error(paste('Impossible:', areas.error[i], 'produces negative area'),name='VennDiagramLogger')
      stop(paste('Impossible:', areas.error[i], 'produces negative area'));
    }
  }
  
  ## rescaling area labels to be proportional to area
  if(length(cex.prop) > 0){
    
    if(length(cex.prop) != 1) {
      flog.error('Value passed to cex.prop is not length 1',name='VennDiagramLogger')
      stop('Value passed to cex.prop is not length 1')
    }
    
    ## figure out what function to use
    func = cex.prop
    if (!is(cex.prop, 'function')) {
      if(cex.prop == 'lin'){
        func = function(x) x
      }
      else if(cex.prop == 'log10'){
        func = log10
      }
      else flog.error(paste0('Unknown value passed to cex.prop: ', cex.prop),name='VennDiagramLogger')
      stop(paste0('Unknown value passed to cex.prop: ', cex.prop))
    }
    
    ## rescale areas
    maxArea = max(areas)            
    for(i in 1:length(areas)){                
      cex[i] = cex[i] * func(areas[i]) / func(maxArea)
      if(cex[i] <= 0) stop(paste0('Error in rescaling of area labels: the label of area ',
                                  i, ' is less than or equal to zero'))
    }
  }
  
  # initialize gList to hold all Grobs generated
  grob.list <- gList();
  
  # plot the ellipses of the Venn diagram
  ellipse.positions <- matrix(
    nrow = 4,
    ncol = 7
  );
  colnames(ellipse.positions) <- c('x', 'y', 'a', 'b', 'rotation', 'fill.mapping', 'line.mapping');
  
  ellipse.positions[1,] <- c(0.65, 0.47, 0.35, 0.20,  45, 2, 2);
  ellipse.positions[2,] <- c(0.35, 0.47, 0.35, 0.20, 135, 1, 1);
  ellipse.positions[3,] <- c(0.50, 0.57, 0.33, 0.15,  45, 4, 4);
  ellipse.positions[4,] <- c(0.50, 0.57, 0.35, 0.15, 135, 3, 3);
  
  # draw the ellipses themselves
  for (i in 1:4) {
    grob.list <- gList(
      grob.list,
      VennDiagram::ellipse(
        x = ellipse.positions[i,'x'],
        y = ellipse.positions[i,'y'],
        a = ellipse.positions[i,'a'],
        b = ellipse.positions[i,'b'],
        rotation = ellipse.positions[i, 'rotation'],
        gp = gpar(
          lty = 0,
          fill = fill[ellipse.positions[i,'fill.mapping']],
          alpha = alpha[ellipse.positions[i,'fill.mapping']]
        )
      )
    );
  }
  
  # draw the ellipse borders
  for (i in 1:4) {
    grob.list <- gList(
      grob.list,
      ellipse(
        x = ellipse.positions[i,'x'],
        y = ellipse.positions[i,'y'],
        a = ellipse.positions[i,'a'],
        b = ellipse.positions[i,'b'],
        rotation = ellipse.positions[i, 'rotation'],
        gp = gpar(
          lwd = lwd[ellipse.positions[i,'line.mapping']],
          lty = lty[ellipse.positions[i,'line.mapping']],
          col = col[ellipse.positions[i,'line.mapping']],
          fill = 'transparent'
        )
      )
    );
  }
  
  # create the labels
  label.matrix <- matrix(
    nrow = 15,
    ncol = 3
  );
  colnames(label.matrix) <- c('label', 'x', 'y');
  
  label.matrix[ 1,] <- c(a1,  0.350, 0.77);
  label.matrix[ 2,] <- c(a2,  0.500, 0.69);
  label.matrix[ 3,] <- c(a3,  0.650, 0.77);
  label.matrix[ 4,] <- c(a4,  0.310, 0.67);
  label.matrix[ 5,] <- c(a5,  0.400, 0.58);
  label.matrix[ 6,] <- c(a6,  0.500, 0.47);
  label.matrix[ 7,] <- c(a7,  0.600, 0.58);
  label.matrix[ 8,] <- c(a8,  0.690, 0.67);
  label.matrix[ 9,] <- c(a9,  0.180, 0.58);
  label.matrix[10,] <- c(a10, 0.320, 0.42);
  label.matrix[11,] <- c(a11, 0.425, 0.38);
  label.matrix[12,] <- c(a12, 0.575, 0.38);
  label.matrix[13,] <- c(a13, 0.680, 0.42);
  label.matrix[14,] <- c(a14, 0.820, 0.58);
  label.matrix[15,] <- c(a15, 0.500, 0.28);
  label.matrix[, 1L] <- signif(x = label.matrix[, 1L] / N_Sets, digits = 4)
  
  processedLabels <- rep('',length(label.matrix[,'label']));
  if(print.mode[1] == 'percent'){
    processedLabels <- paste(signif(label.matrix[,'label']/sum(label.matrix[,'label'])*100,digits=sigdigs),'%',sep='');
    if(isTRUE(print.mode[2] == 'raw'))
    {
      processedLabels <- paste(processedLabels,'\n(',label.matrix[,'label'],')',sep='');
    }
  }
  if(print.mode[1] == 'raw'){
    processedLabels <- label.matrix[,'label'];
    if(isTRUE(print.mode[2] == 'percent'))
    {
      processedLabels <- paste(processedLabels,'\n(',paste(signif(label.matrix[,'label']/sum(label.matrix[,'label'])*100,digits=sigdigs),'%)',sep=''),sep='');
    }
  }
  
  
  for (i in 1:nrow(label.matrix)) {
    grob.list <- gList(
      grob.list,
      textGrob(
        label = processedLabels[i],
        x = label.matrix[i,'x'],
        y = label.matrix[i,'y'],
        gp = gpar(
          col = label.col[i],
          cex = cex[i],
          fontface = fontface[i],
          fontfamily = fontfamily[i]
        )
      )
    );
  }
  
  
  # find the location and plot all the category names
  cat.pos.x <- c(0.18, 0.82, 0.35, 0.65);
  cat.pos.y <- c(0.58, 0.58, 0.77, 0.77);
  
  for (i in 1:4) {
    
    # work out location of the category label
    this.cat.pos <- find.cat.pos(
      x = cat.pos.x[i],
      y = cat.pos.y[i],
      pos = cat.pos[i],
      dist = cat.dist[i]
    );
    
    # then print it
    grob.list <- gList(
      grob.list,
      textGrob(
        label = category[i],
        x = this.cat.pos$x,
        y = this.cat.pos$y,
        just = cat.just[[i]],
        gp = gpar(
          col = cat.col[i],
          cex = cat.cex[i],
          fontface = cat.fontface[i],
          fontfamily = cat.fontfamily[i]
        )
      )
    );
  }
  
  # adjust grob.list to fit and return grob.list
  grob.list <- VennDiagram::adjust.venn(VennDiagram::rotate.venn.degrees(grob.list, rotation.degree, rotation.centre[1], rotation.centre[2]), ...);
  
  # draw diagram before returning gList is specified by user
  if (ind) { grid.draw(grob.list); }
  return(grob.list);
}


GRangeToDFrame <- function(GRangeObject,
                           FeaturesToCollect = c("gene",
                                                 "pseudogene"),
                           Verbose = FALSE) {
  # search until no more searches are necessary
  s1 <- as.character(GRangeObject$type)
  s2 <- GRangeObject$gene
  if (is.null(s2)) {
    s2 <- rep(NA, length(s1))
    s2[s1 %in% FeaturesToCollect] <- paste("unnamed_feature",
                                           seq(sum(s1 %in% FeaturesToCollect)))
  }
  s3 <- GRangeObject@strand
  s4 <- GRangeObject@ranges
  s5 <- as.character(GRangeObject@seqnames)
  s6 <- GRangeObject$Note
  s7 <- GRangeObject$Parent
  s8 <- GRangeObject$ID
  s9 <- !is.na(s2)
  
  # return(list(s1,
  #             s2,
  #             s3,
  #             s4,
  #             s5,
  #             s6,
  #             s7,
  #             s8,
  #             s9))
  
  TOTAL <- sum(table(s1[s1 %in% FeaturesToCollect]))
  
  if (TOTAL == 0) {
    if (Verbose) {
      print("No Features present to collect.")
    }
    return(NULL)
  }
  # print(TOTAL)
  CONTINUE <- TRUE
  KEEP <- vector(mode = "logical",
                 length = length(s1))
  START <- STOP <- vector(mode = "integer",
                          length = length(s1))
  NOTE <- CONTIG <- TYPE <- ID <- NAME <- vector(mode = "character",
                                                 length = length(s1))
  COUNT <- 1L
  FOUNDFEATURES <- 0L
  if (Verbose) {
    pBar <- txtProgressBar(style = 1L)
    TIMESTART <- Sys.time()
  }
  while (CONTINUE) {
    # is the line a line to evaluate
    # check its children
    if (s1[COUNT] %in% FeaturesToCollect) {
      w1 <- which(s7 == s8[COUNT])
      w1 <- which(lengths(w1) > 0L)
      # print(w1)
      # if the feature has any children
      if (length(w1) > 0L) {
        ph1 <- ""
        for (m2 in seq_along(w1)) {
          ph2 <- unlist(s6[w1[m2]])
          # print(ph2)
          if (length(ph2) > 0) {
            if (!is.na(ph2)) {
              # print(nchar(ph1))
              if (nchar(ph1) == 0) {
                ph1 <- ph2
              } else {
                ph1 <- paste(ph1, ph2, sep = "; ")
              }
            }
          } else {
            if (s1[COUNT] == "gene") {
              ph1 <- "normal feature"
            } else {
              ph1 <- "non-coding pseudofeature"
            }
          }
        }
        NOTE[COUNT] <- ph1
      } else {
        # feature has no children, what to do here?
        NOTE[COUNT] <- "child lines absent"
      }
      START[COUNT] <- s4@start[COUNT]
      STOP[COUNT] <- s4@start[COUNT] + s4@width[COUNT] - 1L
      CONTIG[COUNT] <- s5[COUNT]
      TYPE[COUNT] <- s1[COUNT]
      NAME[COUNT] <- s8[COUNT]
      
      if (s9[COUNT]) {
        ID[COUNT] <- s2[COUNT]
      } else {
        ID[COUNT] <- ""
      }
      KEEP[COUNT] <- TRUE
      FOUNDFEATURES <- FOUNDFEATURES + 1L
    }
    
    if (FOUNDFEATURES >= TOTAL) {
      CONTINUE <- FALSE
    } else {
      COUNT <- COUNT + 1L
    }
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = FOUNDFEATURES / TOTAL)
    }
    
  }
  if (Verbose) {
    close(pBar)
    cat("\n")
    TIMEEND <- Sys.time()
    print(TIMEEND - TIMESTART)
  }
  res <- DataFrame("Start" = START[KEEP],
                   "Stop" = STOP[KEEP],
                   "Type" = TYPE[KEEP],
                   "Contig" = CONTIG[KEEP],
                   "ID" = ID[KEEP],
                   "Note" = NOTE[KEEP],
                   "Name" = NAME[KEEP])
  return(res)
}

NoteCheck <- function(NoteVector,
                      CheckVector) {
  Res <- sapply(X = CheckVector,
                FUN = function(y) {
                  sum(grepl(pattern = y,
                            x = NoteVector))
                })
  return(Res)
}

###### -- variables and paths -------------------------------------------------

ASSEMBLERS <- c("Megahit",
                "SKESA",
                "SPADES",
                "Unicycler")

PseudoTypeVector <- c("frameshifted",
                      "internal stop",
                      "partial abbutting assembly gap",
                      "partial in the middle",
                      "missing C-terminus",
                      "missing N-terminus")

ColVec1 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9', '#ffffff', '#000000')
ColVec2 <- paste0(ColVec1,
                  "33")
PATH01 <- "~/Data/20230111_Factorial_Annotations"
PATH02 <- "~/Data/20230111_Factorial_Assemblies"

load(file = "~/Repos/PseudogenePlots/InputData/Neisseria_v02.RData",
     verbose = TRUE)
load(file = "~/Repos/Pseudogenes/SearchResults.RData",
     verbose = TRUE)
load(file = "~/Repos/Pseudogenes/SearchResults2.RData",
     verbose = TRUE)
FILES01 <- list.files(path = PATH01,
                      full.names = TRUE)
FILES02 <- list.files(path = PATH02,
                      full.names = TRUE)
JobMap <- do.call(rbind,
                  list(read.table("~/Repos/20221207_FactorialAssembly/JobMap_Assembly_1.txt"),
                       read.table("~/Repos/20221207_FactorialAssembly/JobMap_Assembly_2.txt"),
                       read.table("~/Repos/20221207_FactorialAssembly/JobMap_Assembly_3.txt")))
o1 <- order(JobMap$V5)
JobMap <- JobMap[o1, ]
colnames(JobMap) <- c("Run",
                      "Trim",
                      "Quality",
                      "Downsample",
                      "PersistentID",
                      "SP_Placeholder")

###### -- code body -----------------------------------------------------------


L01 <- length(ASSEMBLERS)
UBS <- unique(dat1$BioSample)

res01 <- vector(mode = "list",
                length = length(UBS))

for (m1 in seq_along(UBS)) {
  w2 <- which(dat1$BioSample == UBS[m1])
  
  w3 <- which(JobMap$Run == dat1$Run[w2[1]] &
                JobMap$Trim == TRUE &
                JobMap$Quality == 4L &
                JobMap$Downsample == 1L)
  w4 <- JobMap$PersistentID[w3]
  
  FILES03 <- FILES01[grepl(pattern = formatC(x = w4, width = 6, flag = 0, format = "d"),
                           x = FILES01)]
  FILES04 <- FILES02[grepl(pattern = formatC(x = w4, width = 6, flag = 0, format = "d"),
                           x = FILES02)]
  if (length(FILES03) < 4L | 
      length(FILES04) < 4L) {
    print(paste(m1,
                "is missing a required file."))
    next
  }
  seqs01 <- seqs02 <- GC01 <- GC02 <- GC03 <- NOTES <- vector(mode = "list",
                                                              length = length(ASSEMBLERS))
  
  for (m2 in seq_along(ASSEMBLERS)) {
    w6 <- grepl(pattern = ASSEMBLERS[m2],
                x = FILES03)
    w7 <- grepl(pattern = ASSEMBLERS[m2],
                x = FILES04)
    
    seqs01[[m2]] <- readDNAStringSet(FILES04[w7])
    GC01[[m2]] <- rtracklayer::import(FILES03[w6])
    GC02[[m2]] <- gffToDataFrame(GFF = FILES03[w6],
                                 Verbose = TRUE)
    GC01[[m2]] <- GRangeToDFrame(GRangeObject = GC01[[m2]],
                                 Verbose = TRUE)
    
    # ensure names are unique for later
    GC02[[m2]]$ID <- paste0(GC02[[m2]]$ID,
                            "_",
                            m2)
    GC01[[m2]]$Name <- paste0(GC01[[m2]]$Name,
                              "_",
                              m2)
    fsnote <- grepl(x = GC01[[m2]]$Note,
                    pattern = "frameshifted")
    isnote <- grepl(x = GC01[[m2]]$Note,
                    pattern = "internal stop")
    
    # o1 <- match(x = GC01[[m2]]$Name,
    #             table = GC02[[m2]]$ID)
    o3 <- match(x = GC02[[m2]]$ID,
                table = GC01[[m2]]$Name)
    o2 <- unlist(regmatches(x = names(seqs01[[m2]]),
                            m = gregexpr(pattern = "^[^ ]+",
                                         text = names(seqs01[[m2]]))))
    o2 <- match(x = GC02[[m2]]$Contig,
                table = o2)
    GC02[[m2]] <- cbind(GC02[[m2]],
                        "FS" = fsnote[o3],
                        "IS" = isnote[o3],
                        "ContigMax" = width(seqs01[[m2]])[o2])
    NOTES[[m2]] <- NoteCheck(NoteVector = GC01[[m2]]$Note,
                             CheckVector = PseudoTypeVector)
    
    seqs02[[m2]] <- ExtractBy(x = GC02[[m2]][GC02[[m2]]$Coding, ],
                              y = seqs01[[m2]])
  }
  # 
  # x <- c(seqs02[[1]][GC02[[1]]$FS[GC02[[1]]$Coding]],
  #        seqs02[[2]][GC02[[2]]$FS[GC02[[2]]$Coding]],
  #        seqs02[[3]][GC02[[3]]$FS[GC02[[3]]$Coding]],
  #        seqs02[[4]][GC02[[4]]$FS[GC02[[4]]$Coding]])
  
  x <- do.call(c, seqs02)
  y <- Clusterize(myXStringSet = x,
                  cutoff = 0.1,
                  includeTerminalGaps = TRUE,
                  penalizeGapLetterMatches = TRUE)
  z1 <- tapply(X = rownames(y),
               INDEX = y$cluster,
               FUN = c)
  # z1 <- tapply(X = seq_along(y$cluster),
  #              INDEX = y$cluster,
  #              FUN = c)
  # z2 <- lapply(X = z1,
  #              FUN = function(y) {
  #                x[y]
  #              })
  # z3 <- lapply(X = z2,
  #              FUN = function(x) {
  #                if (length(x) >= 2L) {
  #                  AlignSeqs(myXStringSet = x,
  #                            verbose = FALSE)
  #                } else {
  #                  NULL
  #                }
  #              })
  # z4 <- lapply(X = z3,
  #              FUN = function(x) {
  #                if (!is.null(x)) {
  #                  DistanceMatrix(myXStringSet = x,
  #                                 verbose = FALSE,
  #                                 includeTerminalGaps = TRUE)
  #                }
  #              })
  # z5 <- unlist(sapply(X = z4,
  #                     FUN = function(x) {
  #                       if (!is.null(x)) {
  #                         x[upper.tri(x)]
  #                       }
  #                     }))
  
  z5 <- list(names(seqs02[[1]][GC02[[1]]$FS[GC02[[1]]$Coding]]),
             names(seqs02[[2]][GC02[[2]]$FS[GC02[[2]]$Coding]]),
             names(seqs02[[3]][GC02[[3]]$FS[GC02[[3]]$Coding]]),
             names(seqs02[[4]][GC02[[4]]$FS[GC02[[4]]$Coding]]))
  
  z6 <- list(names(seqs02[[1]][GC02[[1]]$IS[GC02[[1]]$Coding]]),
             names(seqs02[[2]][GC02[[2]]$IS[GC02[[2]]$Coding]]),
             names(seqs02[[3]][GC02[[3]]$IS[GC02[[3]]$Coding]]),
             names(seqs02[[4]][GC02[[4]]$IS[GC02[[4]]$Coding]]))
  
  res01[[m1]] <- list("clust" = z1,
                      "FSIDs" = z5,
                      "ISIDs" = z6)
  print(Sys.time())
  print(m1)
}

save(res01,
     file = "~/Data/20230427_NeisseriaGeneSets.RData",
     compress = "xz")

###### -- divy up sets into their respective sets -----------------------------
# a
# a b
# a b c
# a b c d
# b
# b c
# b c d
# c
# c d
# a c
# a c d
# b d
# a d

pBar <- txtProgressBar(style = 1)
PBAR <- length(res01)
res02 <- vector(mode = "list",
                length = PBAR)

for (m1 in seq_along(res01)) {
  if (length(res01[[m1]]) > 0) {
    res02[[m1]] <- matrix(data = F,
                          nrow = length(res01[[m1]][[1]]),
                          ncol = 4L)
    for (m2 in seq_along(res01[[m1]][[1]])) {
      if (any(res01[[m1]][[1]][[m2]] %in% res01[[m1]][[2]][[1]])) {
        res02[[m1]][m2, 1L] <- T
      }
      if (any(res01[[m1]][[1]][[m2]] %in% res01[[m1]][[2]][[2]])) {
        res02[[m1]][m2, 2L] <- T
      }
      if (any(res01[[m1]][[1]][[m2]] %in% res01[[m1]][[2]][[3]])) {
        res02[[m1]][m2, 3L] <- T
      }
      if (any(res01[[m1]][[1]][[m2]] %in% res01[[m1]][[2]][[4]])) {
        res02[[m1]][m2, 4L] <- T
      }
    }
  }
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

pBar <- txtProgressBar(style = 1)
PBAR <- length(res01)
res03 <- vector(mode = "list",
                length = length(res02))
for (m1 in seq_along(res02)) {
  if (length(res02[[m1]]) > 0) {
    res03[[m1]] <- apply(X = res02[[m1]],
                         MARGIN = 2L,
                         FUN = function(x) {
                           paste(m1, which(x), sep = "_")
                         })
  }
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

uids <- unique(unlist(res03))
res03 <- do.call(cbind, res03)
res03 <- lapply(X = 1:4,
                FUN = function(x) {
                  as.integer(factor(x = do.call(c, res03[x, ]),
                                    levels = uids))
                })

res03 <- apply(X = res02[[1]],
               MARGIN = 2L,
               FUN = function(x) {
                 which(x)
               })
names(res03) <- ASSEMBLERS

FSVenn <- res03

pBar <- txtProgressBar(style = 1)
PBAR <- length(res01)
res02 <- vector(mode = "list",
                length = PBAR)

for (m1 in seq_along(res01)) {
  if (length(res01[[m1]]) > 0) {
    res02[[m1]] <- matrix(data = F,
                          nrow = length(res01[[m1]][[1]]),
                          ncol = 4L)
    for (m2 in seq_along(res01[[m1]][[1]])) {
      if (any(res01[[m1]][[1]][[m2]] %in% res01[[m1]][[3]][[1]])) {
        res02[[m1]][m2, 1L] <- T
      }
      if (any(res01[[m1]][[1]][[m2]] %in% res01[[m1]][[3]][[2]])) {
        res02[[m1]][m2, 2L] <- T
      }
      if (any(res01[[m1]][[1]][[m2]] %in% res01[[m1]][[3]][[3]])) {
        res02[[m1]][m2, 3L] <- T
      }
      if (any(res01[[m1]][[1]][[m2]] %in% res01[[m1]][[3]][[4]])) {
        res02[[m1]][m2, 4L] <- T
      }
    }
  }
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

pBar <- txtProgressBar(style = 1)
PBAR <- length(res01)
res03 <- vector(mode = "list",
                length = length(res02))
for (m1 in seq_along(res02)) {
  if (length(res02[[m1]]) > 0) {
    res03[[m1]] <- apply(X = res02[[m1]],
                         MARGIN = 2L,
                         FUN = function(x) {
                           paste(m1, which(x), sep = "_")
                         })
  }
  setTxtProgressBar(pb = pBar,
                    value = m1 / PBAR)
}
close(pBar)

uids <- unique(unlist(res03))
res03 <- do.call(cbind, res03)
res03 <- lapply(X = 1:4,
                FUN = function(x) {
                  as.integer(factor(x = do.call(c, res03[x, ]),
                                    levels = uids))
                })

res03 <- apply(X = res02[[1]],
               MARGIN = 2L,
               FUN = function(x) {
                 which(x)
               })
names(res03) <- ASSEMBLERS

ISVenn <- res03

save(FSVenn,
     ISVenn,
     file = "~/Repos/PseudogenePlots/InputData/Neisseria_v03.RData",
     compress = "xz")

venn.diagram(x = res03,
             filename = "testvenndiagram.png")

venn.diagram.adhoc(x = res03,
                   # N_Sets = sum(lengths(res01) == 3L),
                   sigdigs = 2,
                   N_Sets = 1,
                   print.mode = "percent",
                   filename = "dummystring")
