#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <caml/bigarray.h> /* Added by Lauren */
#include <caml/memory.h> /* Added by Lauren */
#include <caml/custom.h> 
#include <caml/intext.h> 
#include <caml/fail.h> 
#include <caml/alloc.h> 
#include <caml/mlvalues.h> 
#include "invdist.h"
#include "correction.h"
#include "binencode.h"
#include "structs.h"
#include "all_sorting_reversals.h"
#include "lists.h"
#include "med_util.h"
#include "sorting_reversal_median.h"
#include "inversion_median.h"
#include "inversion_median_alberto.h"
#include "growTree.h"
#include "condense.h"
#include "condense3.h"
#include "inittree.h"
#include "labeltree.h"
#include "specialinit.h"
#include "convert.h"
#include "bbtsp.h"

#include "lk_main.h"

#include "mgr.h"

#define Genome_matrix_struct(a) ((struct genome_struct *) Data_custom_val(a))

VertexFactory *newvf = NULL; 
int DOBRANCH;


struct genome_arr_t
{
    struct genome_struct *genome_ptr;
    int num_genome;
    int num_gene;
};

void
ini_mem_4_all (int num_genes )
{
 /*  output_genome =
        ( struct genome_struct * ) calloc ( 1,
                                            sizeof ( struct genome_struct ) );
    output_genome->genes = ( int * ) calloc ( num_genes, sizeof ( int ) );
    output_genome->delimiters = ( int * ) calloc ( num_genes, sizeof ( int ) );
    output_genome->deli_num = 0;
    */
}

void
free_mem_4_all ()
{/*
    free ( output_genome->genes );
    free ( output_genome->delimiters );
    free ( output_genome ); */
}

void grappa_CAML_genome_arr_free (value c_genome_arr) {
    //printf("Start of genone_arr_CAML_free\n"); fflush(stdout); 
    struct genome_arr_t *genome_arr;
    struct genome_struct *genome;
    int i;
    genome_arr = (struct genome_arr_t *) Data_custom_val (c_genome_arr);

    for (i = 0 ; i < genome_arr->num_genome; i++){
       genome = genome_arr->genome_ptr + i;
        if (genome != (struct genome_struct *) NULL) {            
            if(genome->genes != (int*) NULL) 
                free (genome->genes);
            else 
                fprintf(stderr," grappa_CAML_genome_arr_free :\
                        try to free empty genes\n");
            if(genome->gnamePtr != (char*) NULL)
                        free (genome->gnamePtr);
            else
                fprintf(stderr," grappa_CAML_genome_arr_free :\
                        try to free emty gnamePtr\n");
            if(genome->delimiters != (int*) NULL)
                free (genome->delimiters);
            else
                 fprintf(stderr," grappa_CAML_genome_arr_free :\
                        try to free emty delimiters\n");
        }
    }
    free (genome_arr->genome_ptr); 
    //printf("End of genone_arr_CAML_free\n"); fflush(stdout); 
    return;
}

static struct custom_operations genomeOps = {
    "http://www.amnh.org/poy/genome/grappa.0.1",
    custom_finalize_default, 
  //  &genome_CAML_free,
    custom_compare_default,
    custom_hash_default,
    custom_serialize_default,
    custom_deserialize_default,
};

static struct custom_operations genomeArrOps = {
    "http://www.amnh.org/poy/genome/grappa.0.1",
    &grappa_CAML_genome_arr_free, 
//    custom_finalize_default,
    custom_compare_default,
    custom_hash_default,
    custom_serialize_default,
    custom_deserialize_default,
};

value grappa_CAML_print_genome(value c_genome, value c_num_gen) {
    CAMLparam2(c_genome, c_num_gen);
    struct genome_struct *genome;
    int i, num_gen;
    genome = (struct genome_struct *) Data_custom_val(c_genome);
    num_gen = Int_val (c_num_gen);
    printf ("Number gene for printing: %i\n", num_gen); fflush (stdout);
    printf("%s: ",genome->gnamePtr); fflush (stdout);
    for (i = 0 ; i < num_gen ; i++) {
        printf(" %2d ",genome->genes[i]);
    }
    printf("\nEnd of printing the genome\n"); fflush (stdout);
    CAMLreturn (Val_unit);
}
            
value grappa_CAML_print_genome_arr (value c_genome_arr, value c_num_genome, value c_num_gen) {
    CAMLparam3(c_genome_arr, c_num_genome, c_num_gen);
    struct genome_arr_t *genome_arr;
    int i,j;
    int num_genome, num_gen;
    genome_arr = (struct genome_arr_t *) Data_custom_val(c_genome_arr);
    num_genome = Int_val(c_num_genome);
    num_gen = Int_val(c_num_gen);

    printf ("Number genome: %i, number gene: %i\n", num_genome, num_gen);
    for (i = 0; i < num_genome; i++) {
        printf("%s: ",(genome_arr->genome_ptr + i)->gnamePtr);
        printf  ("%p\n", genome_arr->genome_ptr + i); fflush (stdout);
        for (j = 0 ; j < num_gen ; j++)
            printf(" %2d", (genome_arr->genome_ptr + i)->genes[j]);
        printf("\n"); fflush (stdout);
        printf ("\n");
    }
    printf ("\nEnd of all genomes\n"); fflush (stdout);
    CAMLreturn (Val_unit);
}

value grappa_CAML_get_num_genome (value c_genome_arr) {
    int num_genome;
    CAMLparam1 (c_genome_arr);
    struct genome_arr_t *genome_arr;
    genome_arr = (struct genome_arr_t *) Data_custom_val (c_genome_arr);
    num_genome = genome_arr->num_genome;
    CAMLreturn (Val_int(num_genome));
}


value grappa_CAML_get_num_gene (value c_genome_arr) {
    CAMLparam1 (c_genome_arr);
    struct genome_arr_t *genome_arr;
    genome_arr = (struct genome_arr_t *) Data_custom_val (c_genome_arr);
    CAMLreturn (Val_int(genome_arr->num_gene));
}

value grappa_CAML_get_one_genome(value c_genome_arr, value c_index) {
    CAMLparam2(c_genome_arr,c_index);
    struct genome_arr_t *genome_arr;
    struct genome_struct *genome;
    int index;

    CAMLlocal1 (c_genome);

    index = Int_val (c_index);
    genome_arr = (struct genome_arr_t  *) Data_custom_val (c_genome_arr);

    /* Allocate custom block.  This may trigger a garbage collection, causing
     * genome_arr to be moved, so we get the new value of genome_arr. */
    c_genome = alloc_custom(&genomeOps, sizeof (struct genome_struct ), 1, 1000000);
    genome = (struct genome_struct *) Data_custom_val(c_genome);

    genome_arr = (struct genome_arr_t  *) Data_custom_val (c_genome_arr);


/*    int *genes;
    int genome_num;
    char *encoding;
    char *gnamePtr;             
    char parent[32];
*/
    
    genome->genes = (genome_arr->genome_ptr + index)->genes;
    genome->delimiters = (genome_arr->genome_ptr + index)->delimiters;
    genome->deli_num = (genome_arr->genome_ptr + index)->deli_num;
    genome->gnamePtr = (genome_arr->genome_ptr)->gnamePtr;
    genome->encoding = (genome_arr->genome_ptr + index)->encoding;   
    //strcpy(genome->parent, (genome_arr->genome_ptr + index)->parent);
    

    CAMLreturn(c_genome); 
}

value grappa_CAML_get_delimiter_num ( value in_genome )
{
    CAMLparam1(in_genome);
    struct genome_struct * g1;
    g1 = (struct genome_struct *) Data_custom_val (in_genome);
    int delinum;
    delinum = g1->deli_num;
    CAMLreturn(Val_int(delinum)); 
}

value grappa_CAML_get_delimiter_bigarr ( value in_genome, value num_deli )
{
    CAMLparam2(in_genome,num_deli);
    CAMLlocal1(res);
    int NUM_DELI;
    NUM_DELI = Int_val(num_deli);
    struct genome_struct * g1;
    long dims[1]; dims[0] = NUM_DELI;
    g1 = (struct genome_struct *) Data_custom_val (in_genome);
    res =  alloc_bigarray (BIGARRAY_INT32 | BIGARRAY_C_LAYOUT, 1, 
                    g1->delimiters,dims);
    CAMLreturn(res);

}

value grappa_CAML_get_gene_bigarr (value in_genome,value num_gene)
{
    CAMLparam2(in_genome,num_gene);
    CAMLlocal1(res);
    int NUM_GENES;
    NUM_GENES = Int_val(num_gene);
    struct genome_struct * g1;
    long dims[1]; dims[0] = NUM_GENES;
    g1 = (struct genome_struct *) Data_custom_val (in_genome);
    res =  alloc_bigarray (BIGARRAY_INT32 | BIGARRAY_C_LAYOUT, 1, 
                    g1->genes,dims);
    CAMLreturn(res);
}

/* Added by Lauren to take in Ocaml values. 
 Returns the inversion distance between gene1 and gene2.
 gene1-- the first genome
 gene2-- the 2nd genome
 num-- the number of gene in each genome (i.e., the length of gene1;
 the length of both must be the same)
 circular-- an int that's positive if genome are circular */
value grappa_CAML_cmp_inv_dis(value c_gene1, value c_gene2, 
				      value c_num_gen, value c_circular) {
    CAMLparam4(c_gene1, c_gene2, c_num_gen, c_circular); 
    int num_gene, distance, circ, num_chromosome;
    struct genome_struct *g1, *g2;
    g1 = (struct genome_struct *) Data_custom_val (c_gene1);
    g2 = (struct genome_struct *) Data_custom_val (c_gene2);
    num_gene = Int_val (c_num_gen);
    circ = Int_val (c_circular);
    int deli_num1 = g1->deli_num;
    int deli_num2 = g2->deli_num;
  /*  fprintf(stdout,"grappa, c_cmp_inv_dis with deli_num1=%d, deli_num2=%d \n",
            deli_num1, deli_num2);
    fflush(stdout);
    int i;
        fprintf(stdout,"g1 = {");
        for(i=0;i<num_gene;i++)
        { 
            fprintf(stdout,"%d,",g1->genes[i]);
        }
        fprintf(stdout,"} ; g2 = { ");
        for(i=0;i<num_gene;i++)
        { 
            fprintf(stdout,"%d,",g2->genes[i]);
        }fprintf(stdout," }\n");
        fflush(stdout); */
    if ((deli_num1>1)||(deli_num2>1)) //deal with multichromosome input
    {
        distance = mgr_invdist(g1->genes,g2->genes,num_gene,g1->delimiters,g2->delimiters,deli_num1,deli_num2); 
    }
    else
    {
        if (circ == 1) {
             distance = invdist_circular_nomem(g1, g2, num_gene); 
        } else {
            distance = invdist_noncircular_nomem(g1, g2, 0, num_gene);
        }
    }
    CAMLreturn(Val_int(distance));
}


value grappa_CAML_better_capping (value c_gene1, value c_gene2, value num_genes)
{
    CAMLparam3(c_gene1,c_gene2,num_genes);
    int NUM_GENES = Int_val(num_genes);
    long dims[1]; dims[0] = NUM_GENES;
    struct genome_struct *g1, *g2, *out_genome_list;
    g1 = (struct genome_struct *) Data_custom_val (c_gene1);
    g2 = (struct genome_struct *) Data_custom_val (c_gene2);
    out_genome_list =
    (struct genome_struct *) malloc (1*sizeof (struct genome_struct) );
    if ( out_genome_list == ( struct genome_struct * ) NULL )
        failwith ("ERROR: genome_list in grappa_CAML_better_capping is NULL" );
    out_genome_list[0].gnamePtr =
            ( char * ) malloc ( MAX_NAME * sizeof ( char ) );
    sprintf (out_genome_list[0].gnamePtr, "%i", 0);
    if ( out_genome_list[0].gnamePtr == ( char * ) NULL )
            failwith( "ERROR: gname of genome_list in grappa_CAML_better_capping is NULL" );
    out_genome_list[0].genes =( int * ) malloc ( NUM_GENES * sizeof ( int ) );
    out_genome_list[0].delimiters = (int *) malloc (NUM_GENES * sizeof (int) );

    better_capping (g1->genes,g2->genes,NUM_GENES,g1->delimiters,g2->delimiters,g1->deli_num,g2->deli_num,out_genome_list);

    struct genome_arr_t *out_genome_arr;
    CAMLlocal1 (c_genome_arr);
    c_genome_arr = alloc_custom(&genomeArrOps, sizeof(struct genome_arr_t), 1, 10000);
    out_genome_arr = (struct genome_arr_t *) Data_custom_val(c_genome_arr);
    out_genome_arr->genome_ptr = out_genome_list;    
    out_genome_arr->num_genome = 1;
    out_genome_arr->num_gene = NUM_GENES;
    CAMLreturn(c_genome_arr); 

}

value 
grappa_CAML_inv_med 
(value medsov, value c_gene1, value c_gene2, value c_gene3, value num_genes,value circular)
{
    CAMLparam5(medsov,c_gene1,c_gene2,c_gene3,num_genes);
    CAMLxparam1(circular);
    CAMLlocal1(res);
    int MEDIAN_SOLVER;
    struct genome_struct *g1, *g2, *g3;
    struct genome_struct *gen[3];
    struct genome_struct *out_genome_list;
    struct genome_arr_t *out_genome_arr;
    int CIRCULAR;   
    int NUM_GENES;
    int num_cond;
    int old_max_num_genes;
    int multichromosome=0;
    MEDIAN_SOLVER = Int_val(medsov);
    g1 = (struct genome_struct *) Data_custom_val (c_gene1);
    g2 = (struct genome_struct *) Data_custom_val (c_gene2);
    g3 = (struct genome_struct *) Data_custom_val (c_gene3);
    CIRCULAR = Int_val(circular);
    NUM_GENES = Int_val(num_genes);
    long dims[1]; dims[0] = NUM_GENES;

    condense3_mem_t * cond3mem_p; cond3mem_p =  &CONDENSE3_MEM;
    convert_mem_t * convertmem_p; convertmem_p = &CONVERT_MEM;
    old_max_num_genes = cond3mem_p->max_num_genes;


    out_genome_list =
        ( struct genome_struct * ) malloc ( 1 *
                                            sizeof ( struct genome_struct ) );
    if ( out_genome_list == ( struct genome_struct * ) NULL )
        fprintf ( stderr, "ERROR: genome_list NULL\n" );
    out_genome_list[0].gnamePtr =
            ( char * ) malloc ( MAX_NAME * sizeof ( char ) );
    sprintf (out_genome_list[0].gnamePtr, "%i", 0);
    if ( out_genome_list[0].gnamePtr == ( char * ) NULL )
    {
            fprintf ( stderr, "ERROR: gname NULL\n" );
    };
    out_genome_list[0].genes =( int * ) malloc ( NUM_GENES * sizeof ( int ) );
    out_genome_list[0].delimiters = (int *) malloc (NUM_GENES * sizeof (int) );

    if (old_max_num_genes >= NUM_GENES) {}
    else
    {
        free_mem_4_all ();
        ini_mem_4_all (NUM_GENES);
        free_mem_4_invdist (&INVDIST_MEM);
        ini_mem_4_invdist (NUM_GENES);
        free_mem_4_albert ();
        ini_mem_4_albert (NUM_GENES);
        free_mem_4_siepel ();
        ini_mem_4_siepel(NUM_GENES);
        free_mem_4_cond3 ();
        ini_mem_4_cond3 (NUM_GENES);
        free_mem_4_convert();
        ini_mem_4_convert(NUM_GENES);
        free_mem_4_mgr();
        mgr_ini_mem(NUM_GENES); 
        //3 times of original gene size is the worst case for multi-chromosome.
    }
   /* if ((g1->deli_num>0)||(g2->deli_num>0)||(g3->deli_num>0)) 
    {
    //if we are dealing with multichromosome, in the worst case, each loci is a singlechromosome ($a$b$c$....), we will need 3 times of memory-size.
        multichromosome = 1;
        if(old_max_num_genes <3*NUM_GENES)
        {
            fprintf(stdout,"expand memory for multichromosome, old_size = %d < 3*%d\n",old_max_num_genes, NUM_GENES);
            fflush(stdout);
             free_mem_4_mgr();
             mgr_ini_mem(3*NUM_GENES,0);
        }
    }*/

    if(MEDIAN_SOLVER<7)
    {
        condense3 ( g1->genes,
                    g2->genes,
                    g3->genes,
                    cond3mem_p->con_g1->genes,
                    cond3mem_p->con_g2->genes,
                    cond3mem_p->con_g3->genes, 
                    NUM_GENES, &num_cond,
                    cond3mem_p->pred1, cond3mem_p->pred2, 
                    cond3mem_p->picked, cond3mem_p->decode );
        if (num_cond>0)
        {
            gen[0] = cond3mem_p->con_g1;
            gen[1] = cond3mem_p->con_g2;
            gen[2] = cond3mem_p->con_g3;
            switch (MEDIAN_SOLVER)
            {
                case 1: //Alberto Capara median solver
                  if ( CIRCULAR )
                          albert_inversion_median_circular 
                              ( gen,num_cond,cond3mem_p->con_med->genes );
                  else
                          albert_inversion_median_noncircular
                              (gen,num_cond,cond3mem_p->con_med->genes );
                break;
                case 2: //A. Siepel median solver
                  find_reversal_median ( cond3mem_p->con_med, gen, num_cond, &SIEPEL_MEM );
                break;
                case 3: //Exact median solver
                   convert2_to_tsp ( gen[0], gen[1], gen[2], convertmem_p->adjl, convertmem_p->adjp,
                                          num_cond, CIRCULAR );
                   bbtsp ( 2 * num_cond, cond3mem_p->con_med->genes, 
                           FALSE, /* cannot use median that does not exist */
                            gen[0]->genes, gen[1]->genes, gen[2]->genes,
                            convertmem_p->adjl, 
                            convertmem_p->neighbors, 
                            convertmem_p->stack, 
                            convertmem_p->outcycle, 
                            convertmem_p->degree,
                            convertmem_p->otherEnd, 
                            convertmem_p->edges, 
                            CIRCULAR );
                break;
                case 4: //Greedy median solver
                convert2_to_tsp ( gen[0], gen[1], gen[2], convertmem_p->adjl, convertmem_p->adjp,
                                          num_cond, CIRCULAR );
                coalestsp ( 2 * num_cond, cond3mem_p->con_med->genes,FALSE, 
                            gen[0]->genes, gen[1]->genes, gen[2]->genes,
                            convertmem_p->adjl, 
                            convertmem_p->neighbors, 
                            convertmem_p->stack, 
                            convertmem_p->outcycle, 
                            convertmem_p->degree,
                            convertmem_p->otherEnd, 
                            convertmem_p->edges,
                            CIRCULAR );
                break;
                /* case5 and case6 need the CONCORDE package  */
                // http://www.tsp.gatech.edu//concorde/downloads/downloads.htm
#ifdef USE_CONCORDE
                case 5: //SimpleLK TSP median solver 
                     convert_to_tsp ( gen[0], gen[1],
                                     gen[2], num_cond, CIRCULAR,
                                     convertmem_p->weights );
                     greedylk ( 2 * num_cond, convertmem_p->weights, 
                                cond3mem_p->con_med->genes,
                                convertmem_p->incycle, 
                                convertmem_p->outcycle );
                    break;
                case 6: //ChainedLK TSP median solver
                    convert_to_tsp ( gen[0], gen[1],
                                     gen[2], num_cond, CIRCULAR,
                                     convertmem_p->weights );
                    chlinkern ( 2 * num_cond,
                            convertmem_p->weights, 
                            cond3mem_p->con_med->genes,
                            convertmem_p->incycle, convertmem_p->outcycle );
                    break;
#endif
                default:
                    fprintf(stderr, "unknown choice of median solver !\n");
                    break;                
            }
            decode3 ( out_genome_list->genes, cond3mem_p->con_med->genes, 
                      cond3mem_p->pred1, cond3mem_p->decode, num_cond );
        }
        else
        {
            out_genome_list = g1;
        }

    }
    else// MEDIAN_SOLVER == 7, MGR median solver
    {
         mgr_med (g1->genes,g2->genes,g3->genes,g1->delimiters,g2->delimiters,g3->delimiters,g1->deli_num,g2->deli_num,g3->deli_num, NUM_GENES,CIRCULAR,out_genome_list);
  /* debug msg
         fprintf(stdout,"out_genome_list = [");
         int xx=0; 
         for(xx=0;xx<NUM_GENES;xx++)
             fprintf(stdout,"%d,",out_genome_list->genes[xx]);
         fprintf(stdout,"]; delimiters = [");
         for(xx=0;xx<out_genome_list->deli_num;xx++)
             fprintf(stdout,"%d",out_genome_list->delimiters[xx]);
         fprintf(stdout,"]\n");
         fflush(stdout);
        debug msg */
    }
    CAMLlocal1 (c_genome_arr);
    c_genome_arr = alloc_custom(&genomeArrOps, sizeof(struct genome_arr_t), 1, 1000000);
    out_genome_arr = (struct genome_arr_t *) Data_custom_val(c_genome_arr);
    out_genome_arr->genome_ptr = out_genome_list;    
    out_genome_arr->num_genome = 1;
    out_genome_arr->num_gene = NUM_GENES;
    CAMLreturn(c_genome_arr); 

}



value 
grappa_CAML_inv_med_bytecode (value * argv, int argn){
    return  (grappa_CAML_inv_med 
        (argv[0],argv[1], argv[2], argv[3], argv[4], argv[5]));
}



/*****************************************************************/
value grappa_CAML_create_empty_genome_arr(value numgenome, value numgene)
{
    int Numgenome, Numgene, i;
    struct genome_struct *genome_list;
    struct genome_arr_t *genome_arr;
    CAMLparam2(numgenome, numgene);
    
    Numgenome = Int_val(numgenome);
    Numgene = Int_val(numgene);
  /************/

    genome_list =
        ( struct genome_struct * ) malloc ( Numgenome *
                                            sizeof ( struct genome_struct ) );
    if ( genome_list == ( struct genome_struct * ) NULL )
        fprintf ( stderr, "ERROR: genome_list NULL\n" );

    for ( i = 0; i < Numgenome; i++ )
    {
         genome_list[i].gnamePtr =
            ( char * ) malloc ( MAX_NAME * sizeof ( char ) );
        sprintf (genome_list[i].gnamePtr, "%i", i);
        if ( genome_list[i].gnamePtr == ( char * ) NULL )
        {
            fprintf ( stderr, "ERROR: gname NULL\n" );
        };

        genome_list[i].genes =( int * ) malloc ( Numgene * sizeof ( int ) );
        genome_list[i].delimiters =
            (int *) malloc (Numgene * sizeof (int) );
        genome_list[i].deli_num = 0;
    }

    CAMLlocal1 (c_genome_arr);
    c_genome_arr = alloc_custom(&genomeArrOps, sizeof(struct genome_arr_t), 1, 10000);
    genome_arr = (struct genome_arr_t *) Data_custom_val(c_genome_arr);
    genome_arr->genome_ptr = genome_list;    
    genome_arr->num_genome = Numgenome;
    genome_arr->num_gene = Numgene;
    
    CAMLreturn(c_genome_arr); 
}

value grappa_CAML_set (value set_what, value c_genome_arr, value c_genome_no, value c_index, value c_gene_no) {
    struct genome_arr_t *genome_arr;
    int index, genome_no, gene_no, deli_no;
    int set_seq;
    CAMLparam5 (set_what, c_genome_arr, c_genome_no, c_index, c_gene_no);
    genome_arr = (struct genome_arr_t  *) Data_custom_val (c_genome_arr);
    genome_no = Int_val (c_genome_no);
    set_seq = Int_val(set_what);
    index = Int_val (c_index);
    if (1==set_seq) // seq genes
    {
    gene_no = Int_val (c_gene_no);
    (genome_arr->genome_ptr + genome_no)->genes[index] = gene_no;
    }
    else if(0==set_seq) //set delimiters
    {
        deli_no = Int_val (c_gene_no);
        (genome_arr->genome_ptr + genome_no)->delimiters[index] = deli_no;
        (genome_arr->genome_ptr + genome_no)->deli_num ++;
    }
    else { fprintf(stderr,"in grappa_CAML_set, unkown type of set (set_what = 0 or 1, 1 is set sequence, 0 is set delimiters)"); }
    CAMLreturn(Val_unit);
}


/*
 * See the OCaml interface for more information about this function. Notice that
 * it produce an inverted list of transformations, that is corrected in the
 * OCaml side. *)
 */

value 
grappa_CAML_inversions (value genes1, value genes2, 
        value c_num_genes, value dist) {
    CAMLparam4(genes1, genes2, c_num_genes, dist);
    CAMLlocal3(resulttmp, result, r);
    List intermediate_reversals_list;
    int num_genes, i, j, inv_dist;
    struct genome_arr_t *genes1_arr, *genes2_arr;
    struct genome_struct *permutation, *origin;
    int *temp_genes;
    Reversal *rev; Reversal revrev;
    result = Val_int(0); /* We start with the empty list */

    inv_dist = Int_val(dist);

    permutation = (struct genome_struct *) Data_custom_val (genes1);
    origin = (struct genome_struct *) Data_custom_val (genes2);
    /* First one in should be ancestor-- the permutation that you want to
       transform into the descendant (even though "origin" is a confusing
       thing to call descendant) */
    num_genes = Int_val(c_num_genes);

    temp_genes = (int *)malloc(num_genes * sizeof(int));

    if (0 == num_genes) CAMLreturn(result);
    /* Initialize list that will be used to store the sorting reversals
       found between the permutations at each step. */
    init_list(&intermediate_reversals_list, (num_genes + 1) * num_genes,
            sizeof(Reversal *));

    i = 0;

    do {

        clear_list(&intermediate_reversals_list);

        find_all_sorting_reversals(&intermediate_reversals_list, NULL, permutation,
                origin, num_genes, NULL);

        if  (list_size(&intermediate_reversals_list) > 0) {
            revrev = list_get(&intermediate_reversals_list, 0).revelement;
            rev = &revrev;
            copy_with_reversal(temp_genes, permutation->genes, num_genes, rev);
            permcopy(permutation->genes, temp_genes, num_genes);
            r = caml_alloc_tuple(2);
            Store_field(r,0,Val_int(rev->start + 1));
            Store_field(r,1,Val_int(rev->stop));
            resulttmp = caml_alloc(2,0);
            Store_field(resulttmp,0,r);
            Store_field(resulttmp,1,result);
            result = resulttmp;
        }

        i++;

    } while (list_size(&intermediate_reversals_list) > 0);

    fflush(stdout); /* Change so can be stderr, too? */

    CAMLreturn(result);
}

void
grappa_ini_convert_mem (int num_genes)
{
    ini_mem_4_convert(num_genes);
    return;
}


void 
grappa_ini_invdis_mem ( int num_genes)
{
    ini_mem_4_invdist(num_genes);
    return;
}

void 
grappa_ini_cond3_mem ( int num_genes)
{
    ini_mem_4_cond3(num_genes);
    return;
}

void 
grappa_ini_albert_mem (int num_genes)
{
    ini_mem_4_albert (num_genes); 
    return;
}

void 
grappa_ini_siepel_mem (int num_genes)
{
    ini_mem_4_siepel (num_genes); 
    return;
}


void
grappa_ini_mem (int num_genes)
{
    ini_mem_4_all (num_genes);
    return;
}


value 
grappa_CAML_initialize (value max_num_genes) {
    CAMLparam1(max_num_genes);
    grappa_ini_mem (Int_val(max_num_genes));
    grappa_ini_invdis_mem (Int_val(max_num_genes));
    grappa_ini_cond3_mem (Int_val(max_num_genes));
    grappa_ini_albert_mem(Int_val(max_num_genes));
    grappa_ini_siepel_mem (Int_val(max_num_genes));
    grappa_ini_convert_mem(Int_val(max_num_genes));
    mgr_ini_mem(Int_val(max_num_genes));
    CAMLreturn(Val_unit);
}
