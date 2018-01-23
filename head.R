senv <- function(value, var_name){
	assign(deparse(substitute(var_name)), value, envir=globalenv())
}

get_noncanonical_regions <- function(ensembl, canonical_transcripts = NULL){
	if(is.null(canonical_transcripts)){
		return(ensembl)
	}else{
		# canonical transcripts
		ensembl %>% 
			mutate(is_canonical = ensembl_transcript_id %in% canonical_transcripts) ->
			ensembl
		
		# filter genes without canonical transcripts
		ensembl %>% 
			group_by(ensembl_gene_id) %>% 
			mutate(has_canonical=any(is_canonical)) %>% 
			ungroup %>% 
			filter(has_canonical == T) %>% 
			select(-has_canonical) ->
			ensembl
		
		ensembl %>% 
			select(chromosome_name, ensembl_gene_id, ensembl_transcript_id, 
						 genomic_coding = genomic_coding_start, is_canonical) %>% 
			mutate(start_end = 1) -> 
			temp1.1
		
		ensembl %>% 
			select(chromosome_name, ensembl_gene_id, ensembl_transcript_id, 
						 genomic_coding = genomic_coding_end, is_canonical) %>% 
			mutate(start_end = -1) -> 
			temp1.2
		
		temp1.1 %>% 
			bind_rows(temp1.2) ->
			temp2
		
		temp2 %>% 
			filter(is_canonical == F) %>% 
			group_by(chromosome_name, ensembl_gene_id) %>% 
			arrange(genomic_coding) %>% 
			mutate(start_end_cumsum=cumsum(start_end)) %>% 
			mutate(is_coding = start_end_cumsum > 0) %>% 
			mutate(p_k = lag(is_coding, default = F)) %>% 
			mutate(pocz_kon = is_coding & p_k) %>% 
			filter(!pocz_kon) %>% 
			select(chromosome_name, ensembl_gene_id, ensembl_transcript_id, genomic_coding, is_canonical) %>% 
			mutate(start_end = rep(c(1, -1), length.out = n())) %>% 
			arrange(genomic_coding) %>% 
			ungroup() ->
			temp3
		
		temp2 %>% 
			filter(is_canonical == T) %>% 
			bind_rows(temp3) %>% 
			group_by(chromosome_name, ensembl_gene_id) %>% 
			arrange(genomic_coding, desc(is_canonical)) %>% 
			mutate(active_canonical = cumsum(is_canonical * start_end)) %>% 
			mutate(active_noncanonical = cumsum((!is_canonical) * start_end)) %>% 
			distinct(chromosome_name, ensembl_gene_id, genomic_coding, .keep_all=T) %>% 
			mutate(pocz = active_canonical == 0 & active_noncanonical == 1) %>% 
			mutate(kon = (active_canonical == 1 & active_noncanonical == 1) | (active_canonical == 0 & active_noncanonical == 0)) %>% 
			mutate(pocz2 = pocz & lead(kon, default = F)) %>% 
			mutate(kon2 = kon & lag(pocz, default = F)) %>% 
			filter(pocz2 | kon2) %>% 
			ungroup() ->
			temp4
		
		temp4 %>% 
			filter(pocz2) %>% 
			select(chromosome_name, ensembl_gene_id, genomic_coding) %>% 
			group_by(ensembl_gene_id) %>% 
			mutate(nr=row_number()) %>% 
			ungroup() ->
			temp5.1
		
		temp4 %>% 
			filter(kon2) %>% 
			select(chromosome_name, ensembl_gene_id, genomic_coding) %>% 
			group_by(ensembl_gene_id) %>% 
			mutate(nr=row_number()) %>% 
			ungroup() ->
			temp5.2
		
		temp5.1 %>% 
			left_join(temp5.2, by=c("chromosome_name", "ensembl_gene_id", "nr")) %>% 
			rename(genomic_coding_start = genomic_coding.x,
						 genomic_coding_end = genomic_coding.y) %>% 
			mutate(ensembl_gene_id = paste0(ensembl_gene_id, ".", nr)) %>% 
			select(-nr) %>% 
			return()
	}
}