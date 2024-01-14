# The helper functions for sequence-feature-informed simulation


#' SimulateRegions
#'
#' The main function to run interval simulation
#'
#' @param region_real
#' @param n_test
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
SimulateRegions <- function(region_real, n_test,  ... ){

  # pool_simu <- BuildSimuPools(region_real, simu_fuc = SimulatePool, simu_times = 3, cores = 5)
  pool_simu <- BuildSimuPools(region_real, simu_times = n_test, ...)



  ## 2.2 QC
  n_cols <- sapply(pool_simu, NCOL, simplify = T)
  peak_to_keep <- n_cols==max(n_cols)
  valid_region_real <- region_real[peak_to_keep,]
  valid_pool_simu <- pool_simu[peak_to_keep]

  # 2.3 Split pools into peaksets
  region_simu <- SplitSimuPools(valid_pool_simu, simu_times = n_test)

  # 2.4 Output real regions and simulate region
  region_real_simu <-  c(list(region_real), region_simu)

  return(region_real_simu)
}

#' SimulatePool
#'
#' simulation peak pool with dynamic toleraance
#'
#' @param interval Input intervals to be simulated
#' @param permutation How many simulation to try to each sample.
#' Too large will waste the simulation
#' @param n_test How many simulations needed?
#' Need to find a sweet point between n_test and permutation
#' @param tolerence The starting tolerance.
#' The starting tolerance should be high to short sequences such as non-B
# It will probably also be better to use real count when the tolerance is small
#' @param increase_tol_after increase_tol_after trying certain numbers of simulation.
#' smaller number makes the simulation faster
#' @param increase_tol_step
#'
#' @return
#' @export
#'
#' @examples
SimulatePool <- function(interval, permutation = 1e+03, n_test = 100,
                         tolerence = 0.02,
                         increase_tol_after = 5,
                         increase_tol_step = 0.01,
                         within_chr = F
                         ) {
  pg <- interval$pG
  ng <- interval$nG
  match_j = NULL
  simu_count = 0
  genome_chr = hg19.size
  if (within_chr) genome_chr = hg19.size %>% filter(chrom %in% interval$chrom)

  while (is.null(match_j) | NROW(match_j) < n_test) {
    j <- valr::bed_random(genome_chr,
                          interval$length-1,
                          n = permutation
                          ) %>%
      CalculateG() %>% arrange(desc(simu_ng)) %>%
      mutate(simu_pg = simu_ng/length) %>%
      mutate(real_ng = ng)

    simu_count = simu_count + 1
    if (simu_count != 1 & simu_count %% increase_tol_after == 0) {
      tolerence = tolerence + increase_tol_step
      message(sprintf("Tolerence increase to %f", tolerence))
    }
    add_j <- j %>%
      filter(between(simu_pg, pg - tolerence, pg + tolerence)) %>%
      mutate(tolerence = tolerence)
    match_j <- match_j %>% bind_rows(add_j) %>% distinct()

    message(sprintf("Simulated %i intervals ...", nrow(match_j)))
  }

  if("peak_id" %in% colnames(interval)) {
    match_j <- mutate(match_j,  peak_id = interval$peak_id)
  }

  return(match_j)
}


#' SimulatePoolCounts
#'
#' # for shorter sequencne
#' # when the seqeunces are short (such as non-B),
#' # instead of using by pap
#'
#'
#' @param interval
#' @param permutation
#' @param n_test
#' @param tolerence
#' @param increase_tol_after
#' @param increase_tol_step
#'
#' @return
#' @export
#'
#' @examples
#'
#' SimulatePoolCounts(region_real[1,], n_test = 100,tolerence = 0, increase_tol_after =1)

SimulatePoolCounts <- function(interval, permutation = 1e+03, n_test = 100,
                         tolerence = 1,
                         increase_tol_after = 5,
                         increase_tol_step = 1
) {
  pg <- interval$pG
  ng <- interval$nG
  match_j = NULL
  simu_count = 0

  while (is.null(match_j) | NROW(match_j) < n_test) {
    j <- valr::bed_random(hg19.size, interval$length-1,n = permutation) %>%
      CalculateG() %>% arrange(desc(simu_ng)) %>%
      mutate(real_ng = ng)


    # increase tolerance
    simu_count = simu_count + 1
    if (simu_count != 1 & simu_count %% increase_tol_after == 0) {
      tolerence = tolerence + increase_tol_step
      message(sprintf("Tolerence increase to %.0f bp", tolerence))
    }

    # fitler and append simulation
    add_j <- j %>%
      filter(between(simu_ng, ng - tolerence, ng + tolerence)) %>%
      mutate(tolerence = tolerence)
    match_j <- match_j %>% bind_rows(add_j) %>% distinct()
    message(sprintf("Simulated %i intervals ...", nrow(match_j)))



  }

  if("peak_id" %in% colnames(interval)) {
    match_j <- mutate(match_j,  peak_id = interval$peak_id)
  }

  return(match_j)
}





#' SimulatePoolCHR
#'
#' simulation peak pool with dynamic toleraance by Chromosome
#'
#' @param interval
#' @param batch_size
#' @param n_test
#' @param tolerence
#' @param increase_tol_after
#'
#' @return
#' @export
#'
#' @examples
SimulatePoolCHR <- function(interval, batch_size = 1e+03, n_test = 100,
                            tolerence = 0.02, increase_tol_after = 5) {
  pg <- interval$pG
  ng <- interval$nG
  match_j = NULL
  simu_count = 0
  interval_reps <- purrr::map_dfr(seq_len(batch_size), ~interval)


  while (is.null(match_j) | NROW(match_j) < n_test) {
    j <- bed_shuffle(interval_reps, hg19.size, within = T) %>%
      CalculateG() %>% arrange(desc(simu_ng)) %>%
      mutate(simu_pg = simu_ng/length) %>%
      mutate(real_ng = ng)

    simu_count = simu_count + 1
    if (simu_count %% increase_tol_after == 0) {
      tolerence = tolerence + 0.01
      message(sprintf("Tolerence increase to %f", tolerence))
    }
    add_j <- j %>%
      filter(between(simu_pg, pg - tolerence, pg + tolerence)) %>%
      mutate(tolerence = tolerence)
    match_j <- match_j %>% bind_rows(add_j) %>% distinct()

    message(sprintf("Simulated %i intervals ...", nrow(match_j)))
  }

  if("peak_id" %in% colnames(interval)) {
    match_j <- mutate(match_j,  peak_id = interval$peak_id)
  }

  return(match_j)
}

#' BuildSimuPools
#'
#' Build Simumation Pool using differecnt simumation functions
#'
#' @param region_real
#' @param simu_fuc  SimulatePool or SimulatePoolCHR
#' @param simu_times
#' @param cores
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
BuildSimuPools <- function(region_real,
                           simu_fuc = SimulatePool,
                           simu_times = 3,
                           cores = 5,
                           ...){


  region_simu <- parallel::mclapply(
    1:nrow(region_real), mc.cores = cores,
    function(i){
      output <- simu_fuc(interval = region_real[i,], n_test = simu_times, ...)
      message(sprintf("Finished (%i)", i))
      return(output)
    }
  )

  if("peak_id" %in% colnames(region_real)) {
    names(region_simu) = region_real$peak_id
  }

  return(region_simu)
}
#' SplitSimuPools
#'
#' This function extracts require simulated genomic region for the pool
#'
#' @param pool_simu  a pool of intervals
#' @param simu_times number of simulation neeed
#'
#' @return
#' @export
#'
#' @examples
SplitSimuPools <- function(pool_simu, simu_times = 100){

  if ("list" %in% class(pool_simu)) {
    pools_simu_table <- pool_simu %>% bind_rows()
    if (! "peak_id" %in% colnames(pools_simu_table)){
      pools_simu_table <- pool_simu %>% bind_rows(.id = "peak_id")
    }
  } else if ("data.frame" %in% class(pool_simu)) {
    pools_simu_table <- pool_simu
  }


  regions_simu <-
    pools_simu_table %>%
    group_by(peak_id) %>%
    dplyr::slice(1:simu_times) %>%
    mutate(index = 1:simu_times) %>%
    group_by(index) %>% group_split()

  return(regions_simu)
}

#' CalculateG
#'
#' @param ranges
#'
#' @return
#' @export
#'
#' @examples
CalculateG <- function(ranges) {
  Granges<- GenomicRanges::makeGRangesFromDataFrame(data.frame(ranges))
  #Get the sequences and compute the GC content
  freqs <- Biostrings::alphabetFrequency(
    getSeq(BSgenome.Hsapiens.UCSC.hg19, Granges))
  output <- ranges %>%
    mutate(simu_ng = freqs[,'G'], length = rowSums(freqs))
  return(output)
}



#' SimuPosistion
#'
#' @param real
#' @param genome
#' @param n_simu
#'
#' @return
#' @export
#'
#' @examples
SimuPosistion <- function(real, genome = hg19.size, n_simu = 5, within = F, add_ng = T){

  simu =
    map(1:n_simu,
        ~ valr::bed_shuffle(x = real,
                            genome = genome,
                            within = within,
                            seed = .x)
        )

  if (add_ng) {
    simu_with_nG =
      simu %>%
      map(~CalculateG(.x) %>% mutate(simu_pg = simu_ng / length))
    return(simu_with_nG)
  }

  return(simu)
}





