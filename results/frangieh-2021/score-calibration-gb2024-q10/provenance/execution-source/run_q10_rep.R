#!/usr/bin/env Rscript
## Extend the frozen GB2024 Frangieh benchmark to adaptive BH q=.10.
## Pair definitions, null fits, reference draws and source statistics are preserved.
for (name in c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")) do.call(Sys.setenv,setNames(list("1"),name))
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
here <- dirname(normalizePath(sub("^--file=","",arg)))
repo <- normalizePath(file.path(here,"../../.."))
source(file.path(repo,"reproduction-code/reproduction-helpers.R"))
.libPaths(unique(c(reproduction_library_paths(file.path(repo,"reproduction-code")),.libPaths())))
suppressPackageStartupMessages({library(pcmt);library(digest);library(jsonlite)})
stopifnot("adaptive_fdr_bh"%in%getNamespaceExports("pcmt"))
source(file.path(repo,"reproduction-code/real-data-analysis/method_config.R"))
source(file.path(repo,"reproduction-code/real-data-analysis/score_helpers.R"))
args<-commandArgs(TRUE);stopifnot(length(args)>=1L)
rep_idx<-as.integer(args[[1L]]);stopifnot(!is.na(rep_idx),rep_idx%in%1:20)
opt<-function(flag,default){j<-match(flag,args);if(is.na(j))default else {stopifnot(j<length(args));args[[j+1L]]}}
cache_dir<-normalizePath(opt("--cache-dir",file.path(repo,"results/frangieh-q10")))
out<-opt("--results-dir",cache_dir);dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out)
pilot_b<-as.integer(opt("--pilot-b-max","0"));stopifnot(pilot_b==0L || pilot_b>=15L)
B_MAX<-if(pilot_b>0L)pilot_b else 25000L;H_STOP<-15L;Q<-.10
sha<-function(path)digest::digest(file=path,algo="sha256",serialize=FALSE)
record<-function(path)list(path=normalizePath(path),bytes=unname(file.info(path)$size),sha256=sha(path))
atomic_json<-function(object,path){tmp<-paste0(path,".tmp");write_json(object,tmp,pretty=TRUE,auto_unbox=TRUE,na="null",null="null",digits=NA);stopifnot(file.rename(tmp,path))}
atomic_rds<-function(object,path){tmp<-paste0(path,".tmp");saveRDS(object,tmp);stopifnot(file.rename(tmp,path))}
verify_records<-function(rows)for(r in rows)stopifnot(file.exists(r$path),file.info(r$path)$size==r$bytes,sha(r$path)==r$sha256)
baseline_dir<-file.path(repo,"reproduction-code/results/frangieh-2021/score-calibration-gb2024-single-guide-v3")
baseline_path<-file.path(baseline_dir,sprintf("score_fdr_rep%03d.rds",rep_idx))
paths<-c(file.path(cache_dir,"cache.rds"),file.path(cache_dir,"cache-manifest.json"),baseline_path)
stopifnot(all(file.exists(paths)))
inputs<-lapply(paths,record)
sources<-lapply(c(file.path(here,"run_q10_rep.R"),file.path(repo,"reproduction-code/real-data-analysis",c("score_helpers.R","method_config.R"))),record)
config<-list(schema_version="frangieh_score5_q10_extension_v1",rep_idx=rep_idx,q=.10,h=H_STOP,B_max=B_MAX,pilot=pilot_b>0L,
  seed_reference_q=.05,seed_reference_q_levels=c(.01,.05),seed_reference_B_max=25000L,
  seed_scope="Exact original q=.05 coordinate streams; q=.10 changes adaptive decisions only. Scalar paired-reference memoization skips redundant deterministic computations.",
  RNGkind=c("Mersenne-Twister","Inversion","Rejection"),pcmt=as.character(packageVersion("pcmt")),
  pcmt_path=find.package("pcmt"),adaptive_body_sha256=digest(body(pcmt::adaptive_fdr_bh),algo="sha256"),
  method_levels=PERTURB_METHOD_LEVELS,n_positive=50L,n_null=950L,n_total=1000L)
fingerprint<-digest(list(config=config,inputs=inputs,sources=sources),algo="sha256")
path<-file.path(out,sprintf("score_fdr_q10_rep%03d.rds",rep_idx));receipt_path<-file.path(out,sprintf("rep%03d-receipt.json",rep_idx))
if(file.exists(path)||file.exists(receipt_path)){
  stopifnot(file.exists(path),file.exists(receipt_path));old<-read_json(receipt_path,simplifyVector=FALSE)
  stopifnot(old$status=="complete",identical(old$fingerprint,fingerprint));verify_records(old$outputs)
  cat("Verified completed q=.10 replicate ",rep_idx," without rewriting outputs.\n",sep="");quit(save="no",status=0)
}
lock<-paste0(path,".lock");stopifnot(dir.create(lock,showWarnings=FALSE))
progress_path<-file.path(out,sprintf("progress-rep%03d.json",rep_idx))
began<-proc.time()[[3L]]
progress<-function(phase,...){atomic_json(c(list(status="running",rep=rep_idx,phase=phase,pid=Sys.getpid(),updated_utc=format(Sys.time(),tz="UTC",usetz=TRUE)),list(...)),progress_path)}
compact<-function(x)x[c("p_value","rejected","futility","active","K","stop_round","t_total","R_threshold")]
fallback<-function(res,idx,p){if(!length(idx))return(res);res$p_value[idx]<-p[idx];res$active<-setdiff(res$active,idx);res$futility<-setdiff(res$futility,idx);r<-idx[p[idx]<=res$R_threshold];res$rejected<-sort(unique(c(setdiff(res$rejected,idx),r)));res}

tryCatch({
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  progress("loading_cache");cache<-readRDS(paths[[1L]]);base<-readRDS(baseline_path)
  pair_table<-base$result$pairs;stopifnot(nrow(pair_table)==1000L,sum(pair_table$is_alt)==50L)
  cp<-cache$canonical_pairs[[sprintf("rep%03d",rep_idx)]]
  stopifnot(identical(pair_table,cp))
  pair_caches<-vector("list",1000L)
  for(j in seq_len(1000L)){
    row<-pair_table[j,];g<-row$gene
    if(isTRUE(row$is_alt)){
      x<-as.integer(cache$pos_X_full[[g]][cache$pos_cell_masks[[g]]]);key<-paste0("pos:",g,"::",g)
    }else{
      ids<-strsplit(row$nt_subset,",",fixed=TRUE)[[1L]]
      x<-as.integer(rowSums(cache$nt_mat_full[cache$pure_nt_idx,ids,drop=FALSE])>=1L);key<-paste0("null::",g)
    }
    gc<-cache$gene_caches[[key]];stopifnot(!is.null(gc),length(x)==row$n_cells,sum(x)==row$n_exposed)
    pc<-precompute_score_pair(x,gc)
    for(k in c("score_sum","information_denom_sq","empirical_denom_sq","information_stat","empirical_stat","score_eff_n","max_score_share"))stopifnot(abs(pc[[k]]-row[[k]])<=1e-7)
    pair_caches[[j]]<-list(X=x,n_cells=length(x),gene_cache=gc,pair_score=pc)
  }
  progress("observed_scores_verified",pairs=1000L)
  asym_idx<-which(vapply(pair_caches,function(v)isTRUE(v$gene_cache$degenerate),logical(1L)))
  information<-pair_table$information_stat;empirical<-pair_table$empirical_stat;raw<-pair_table$score_sum
  observed<-list(parametric_score=information,information_studentized_permutation=information,studentized_permutation=empirical,hgf_signflip_mc=raw,standardized_signflip=raw)
  parametric_p<-base$result$methods$parametric_score$p_value
  stopifnot(max(abs(parametric_p-2*pnorm(abs(information),lower.tail=FALSE)))<1e-12)
  methods<-lapply(observed,function(v)list(observed=v))
  methods$parametric_score$p_value<-parametric_p;methods$parametric_score$rejected<-list(`0.1`=which(p.adjust(parametric_p,"BH")<=Q))
  timing<-list();reference_counts<-list(permutation_computed=0L,permutation_reused=0L,flip_computed=0L,flip_reused=0L)
  draw_seed<-function(t,j,stream)perturb_resample_seed(rep_idx,20L,.05,c(.01,.05),t,j,1000L,25000L,stream)
  ## Each pair has its own hash table of scalar references. Never cache by gene.
  for(stream in c("permutation","flip")){
    memo<-lapply(seq_len(1000L),function(j)new.env(hash=TRUE,parent=emptyenv()))
    local_methods<-if(stream=="permutation")c("information_studentized_permutation","studentized_permutation")else c("hgf_signflip_mc","standardized_signflip")
    for(method in local_methods){
      cat(sprintf("Rep%03d %s q=.10 h=15 Bmax=%d\n",rep_idx,method,B_MAX));flush.console();start<-proc.time()[[3L]]
      which_stat<-if(method%in%c("information_studentized_permutation","hgf_signflip_mc"))1L else 2L
      callback<-function(active_idx,t){
        if(t==1L||t%%500L==0L)progress(method,round=t,active_pairs=length(active_idx),q=Q)
        key<-as.character(t)
        vapply(active_idx,function(j){
          env<-memo[[j]]
          if(exists(key,envir=env,inherits=FALSE)){
            label<-paste0(stream,"_reused");reference_counts[[label]]<<-reference_counts[[label]]+1L
            return(get(key,envir=env,inherits=FALSE)[[which_stat]])
          }
          p<-pair_caches[[j]];set.seed(draw_seed(t,j,if(stream=="permutation")1L else 2L))
          if(stream=="permutation"){
            xp<-p$X[sample.int(p$n_cells)];components<-efficient_score_components(xp,p$gene_cache)
            stats<-c(components$information_stat,components$empirical_stat)
          }else{
            signs<-rademacher_signs(p$n_cells)
            stats<-c(hgf_signflip_stat(signs,p$pair_score),standardized_signflip_stat(signs,p$pair_score))
          }
          assign(key,stats,envir=env);label<-paste0(stream,"_computed");reference_counts[[label]]<<-reference_counts[[label]]+1L
          stats[[which_stat]]
        },numeric(1L))
      }
      adaptive<-pcmt::adaptive_fdr_bh(T_orig=observed[[method]],permute_stat=callback,alpha=Q,h=H_STOP,B_max=B_MAX,alternative="two.sided")
      adaptive<-fallback(adaptive,asym_idx,parametric_p)
      stopifnot(length(adaptive$p_value)==1000L,all(is.finite(adaptive$p_value)),all(adaptive$p_value>=0&adaptive$p_value<=1),!anyDuplicated(adaptive$rejected))
      methods[[method]]$by_q<-list(`0.1`=compact(adaptive));timing[[method]]<-proc.time()[[3L]]-start
      cat(sprintf("Done %s: %.1fs rounds=%d TP=%d FP=%d active=%d\n",method,timing[[method]],adaptive$t_total,sum(pair_table$is_alt[adaptive$rejected]),sum(!pair_table$is_alt[adaptive$rejected]),length(adaptive$active)));flush.console()
    }
    rm(memo);invisible(gc())
  }
  verify_records(inputs);verify_records(sources)
  result<-list(rep=rep_idx,seed=base$result$seed,pairs=pair_table,methods=methods,timings=timing)
  output<-list(config=config,result=result,provenance=list(fingerprint=fingerprint,inputs=inputs,sources=sources,baseline_schema=base$config$schema_version,baseline_provenance=base$provenance,reference_counts=reference_counts,elapsed_seconds=proc.time()[[3L]]-began,session_info=capture.output(sessionInfo())))
  atomic_rds(output,path)
  atomic_json(list(status="complete",rep=rep_idx,fingerprint=fingerprint,config=config,inputs=inputs,sources=sources,outputs=list(record(path)),reference_counts=reference_counts,elapsed_seconds=proc.time()[[3L]]-began,completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE)),receipt_path)
  atomic_json(list(status="complete",rep=rep_idx,updated_utc=format(Sys.time(),tz="UTC",usetz=TRUE)),progress_path)
  unlink(lock,recursive=TRUE);cat(sprintf("Rep%03d q=.10 complete in %.1fs\n",rep_idx,proc.time()[[3L]]-began))
},error=function(e){atomic_json(list(status="failed",rep=rep_idx,error=conditionMessage(e),updated_utc=format(Sys.time(),tz="UTC",usetz=TRUE)),progress_path);unlink(lock,recursive=TRUE);stop(e)})
