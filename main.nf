#!/usr/bin/env nextflow
nextflow.enable.dsl = 2



params.dummy_file = "/home/nehleh/assets/NO_FILE"
dummy_file = file(params.dummy_file)




frequencies = Channel.value('0.2184,0.2606,0.3265,0.1946' )
rates =  Channel.value('0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ,1')
iteration = Channel.value(1..5)
nu_sim = Channel.of(0.05) //0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1
recomlen = Channel.of(100) //100,200,300,400,500,1000,2000,3000,4000,5000
recomrate = Channel.of(0.01) //0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1


params.genome = 10
params.genomelen = '5000'
params.tMRCA = '0.01'
params.analyse = false
params.sim_stat = 1 //0 is just leaves, 1 is for both internal nodes and leaves and 2 is just internal nodes
params.threshold = 0.9
params.seq = "/shared/homes/13298409/SE_2018-20_outbreak_ten.fst"




def helpMessage() {
  log.info"""

  Usage:

  The typical command for running the pipeline is as follows:
    nextflow run main.nf --mode [sim/emp] --other_options

  Options specifying the evolutionary parameters for simulation:
      --mode sim
    Optional arguments:
      --genome               value > 0 (default 10)                 Number of simulated genomes
      --genomelen            value > 1000 (default 100000)          Length of simulated alignment
      --recomlen             value >0 (default 500)                 Length of recombination
      --tMRCA                tMRCA >0 (default 0.01)                TMRCA
      --recomrate            value >0 (default 0.05)                Recombination rate
      --nu_sim               value >0 (default 0.05)                nu value using in simulated data
      --method               pb,cfml,gub (default pb)               Recombination detection methods(PhiloBacteria,ClonalFrameML,Gubbins)
      --analyse              true or false
  Options for detecting recombination in empirical sequence alignments:
    Mandatory arguments:
      --mode emp
      --seq                  fasta file                             Path to input .fasta file
      --method               pb,cfml,gub (default pb)               Recombination detection methods(PhiloBacteria,ClonalFrameML,Gubbins)

  Output Options:
      --outDir               directory                              Output directory to place final output
  --help                                                            This usage statement

   """.stripIndent()
 }



process MakeClonalTree {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_$filename" }
     maxForks 1

     input:
         each iteration

     output:
         tuple val(iteration), path('Clonaltree.tree'), emit: Clonaltree
     """
       make_clonaltree.py -n ${params.genome}  -t ${params.tMRCA}
     """
}

process BaciSim {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1

     input:
         tuple val(iteration), path('Clonaltree')
         each nu_sim
         each recomlen
         each recomrate


     output:
         tuple val(iteration) , path("BaciSimTrees.tree") , emit: BaciSimtrees
         tuple val(iteration) , path('BaciSim_Log.txt') , emit: Recomlog
         tuple val(iteration) , path('BaciSim_Recombination.jpeg'), emit: SimFig
         tuple val(iteration) , path('Sim_nu.csv') , emit: Sim_nu , optional: true
         path('Clonaltree'), emit: Clonaltree
         val iteration , emit: Range
         tuple val(iteration) , path("Recom_stat.csv") , emit: Recomstat
         val nu_sim , emit : nu_sim
         val recomlen , emit : recomlen
         val recomrate , emit: recomrate
         tuple val(iteration) ,path('unroot_Clonaltree.tree'), emit: unroot_Clonaltree

     """
       BaciSim.py -cl ${Clonaltree} -n ${params.genome} -g ${params.genomelen} -l ${recomlen} -r ${recomrate}  -nu ${nu_sim} -s ${params.sim_stat}
     """
}


process Seq_gen {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1
    //errorStrategy 'ignore'

    input:

        tuple val(iteration), path('BaciSimTrees.tree')
        val f
        val r
        val nu_sim
        val recomlen
        val recomrate

    output:
        path "Wholegenome.fasta" , emit: Wholegenome

    """
     numTrees=\$(wc -l < BaciSimTrees.tree | awk '{ print \$1 }')
     seq-gen  -mGTR  -l${params.genomelen} -r$r -f$f -s1 -of BaciSimTrees.tree -p \$numTrees > Wholegenome.fasta
    """
}

process Get_raxml_tree {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1
    errorStrategy 'ignore'

    input:
        path Wholegenome
        val iteration
        val nu_sim
        val recomlen
        val recomrate

    output:
        path 'RAxML_bestTree.tree', emit: MyRaxML
    """
    raxmlHPC -m GTRGAMMA   -p 12345 -s ${Wholegenome} -N 10 -n tree
    """
}

process PhiloBacteria {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1

     input:
        path Clonaltree
        tuple val(iteration), path('Recomlog')
        path Wholegenome
        path MyRaxML
        val nu_sim
        val recomlen
        val recomrate
        tuple val(iteration), path ('Recomstat')

     output:
        path 'Recom_prob_two.h5'                , emit: recom_prob_two   , optional: true
        path 'PB_nu_two.txt'                    , emit: PB_nu_two        , optional: true
        path 'RMSE_PB_two.csv'                  , emit :PB_RMSE_two      , optional: true
        path 'PB_Recom_two.jpeg'                , emit: PB_Recom_two     , optional: true
        path 'PB_Log_two.txt'                   , emit: PB_Log_two       , optional: true
        path 'PB_rcount_two.csv'                , emit: PB_rcount_two    , optional: true
        path 'baci_rcount.csv'                  , emit: Baci_rcount      , optional: true
        path 'baci_delta.csv'                   , emit: Baci_Delta       , optional: true
        path 'PB_delta_two.csv'                 , emit: PB_Delta_two     , optional: true
        path 'PhiloBacter.tree'                 , emit: PB_tree          , optional: true
        path 'PB_dist.csv'                      , emit: PB_dist          , optional: true
        path 'PB_WRF_distance.csv'              , emit: PB_wrf           , optional: true

     """
       phyloHmm.py -t ${MyRaxML}  -a ${Wholegenome}  -cl ${Clonaltree} -rl ${Recomlog} -sim ${params.simulation}  -rs ${Recomstat}

     """
}



process CFML {
   publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
   maxForks 1

     input:
        path Wholegenome
        path MyRaxML
        val iteration
        val nu_sim
        val recomlen
        val recomrate
    output:
        path "CFML.labelled_tree.newick"    , emit: CFMLtree
        path "CFML.importation_status.txt"  , emit: CFML_recom
        path "CFML.result.txt"              , emit: CFML_result
        path "CFML.em.txt"                  , emit: CFML_em

    """
     ClonalFrameML ${MyRaxML} ${Wholegenome} CFML >  CFML.result.txt
    """
}



process CFML_result {
     publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
     maxForks 1

     input:
        path Wholegenome
        path CFML_recom
        path CFMLtree
        tuple val(iteration), path('Recomlog')
        path Clonaltree
        val iteration
        val nu_sim
        val recomlen
        val recomrate

     output:
        path 'CFML_Recombination.jpeg' , emit: CFMLFig
        path 'RMSE_CFML.csv'           , emit : RMSE_CFML  ,optional: true
        path 'CFML_rcount.csv'         , emit: CFML_rcount ,optional: true
        path 'CFML_delta.csv'          , emit: CFML_delta  ,optional: true
        path 'CFML_dist.csv'           , emit: CFML_dist   ,optional: true
        path 'CFML_WRF_distance.csv'   , emit: CFML_wrf    ,optional: true

     """
       CFML_result.py  -cl ${Clonaltree}  -a ${Wholegenome} -cfl ${CFML_recom}  -cft ${CFMLtree}  -rl ${Recomlog} -sim ${params.simulation}
     """
}


process Run_Gubbins {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1

     input:
        path Wholegenome
        val iteration
        val nu_sim
        val recomlen
        val recomrate

    output:
        path "gubbins.node_labelled.final_tree.tre" , emit: Gubbinstree
        path "gubbins.recombination_predictions.gff" , emit: GubbinsRecom
        path "gubbins.per_branch_statistics.csv" , emit: GubbinsStat
    """
     run_gubbins.py  --model GTRGAMMA  --first-model GTRGAMMA -t raxml  -p gubbins   ${Wholegenome}
    """
}

process Gubbins_result {
    publishDir "${params.outDir}" , mode: 'copy' , saveAs:{ filename -> "num_${iteration}/num_${iteration}_nu_${nu_sim}_Rlen_${recomlen}_Rrate_${recomrate}_$filename" }
    maxForks 1
    errorStrategy 'ignore'
     input:
        //tuple val(iteration), path('unroot_Clonaltree')
        path Clonaltree
        tuple val(iteration), path('Recomlog')
        path Wholegenome
        path Gubbinstree
        path GubbinsRecom
        path GubbinsStat
        val nu_sim
        val recomlen
        val recomrate

    output:
        path 'Gubbins_Recombination.jpeg'    , emit: GubbinsFig     , optional: true
        path 'RMSE_Gubbins.csv'              , emit: Rmse_Gubbins   , optional: true
        path 'Gubbinstree_rescale.tree'      , emit: GubbinsRescaletree
        path 'Gubb_rcount.csv'               , emit: Gubb_rcount    , optional: true
        path 'Gubbins_delta.csv'             , emit: Gubb_delta     , optional: true
        path 'Gubb_dist.csv'                 , emit: Gubb_dist      , optional: true
        path 'Gubb_WRF_distance.csv'         , emit: Gubb_wrf       ,optional: true


    """
     Gubbins_result.py  -cl ${Clonaltree} -a ${Wholegenome}  -rl ${Recomlog}  -gl ${GubbinsRecom} -gt ${Gubbinstree} -gs ${GubbinsStat} -sim ${params.simulation}
    """
}



process RMSE_summary {

     publishDir "${PWD}/Summary_${params.outDir}", mode: "copy"
     maxForks 1

     input:
      path  CollectedRMSE_CFML
      path  CollectedRMSE_Gubb
      path  CollectedRMSE_PB_two

     output:
         path 'RMSE_summary.jpeg' , emit: rmse_plot

     """
       rmse_plot_states.py -t rmse_PB_two.csv  -g rmse_Gubbins.csv -c rmse_CFML.csv
     """
}


process RecomCount {

     publishDir "${PWD}/Summary_${params.outDir}", mode: "copy"
     maxForks 1

     input:
      path  CollectedRcount_Baci
      path  CollectedRcount_PB_two
      path  CollectedRcount_Gubb
      path  CollectedRcount_CFML

     output:
         path 'rcount_scatter.jpeg' , emit: rcount_scatter , optional: true
         path 'rcount_boxplot.jpeg' , emit: rcount_boxplot , optional: true

     """
       rcount_plot.py  -t rcount_PB_two.csv  -b rcount_baci.csv  -g rcount_Gubbins.csv  -c rcount_CFML.csv
     """
}

process Delta_Summary {

     publishDir "${PWD}/Summary_${params.outDir}", mode: "copy"
     maxForks 1

     input:
      path  CollectedDelta_Baci
      path  CollectedDelta_PB_two
      path  CollectedDelta_Gubb
      path  CollectedDelta_CFML

     output:
         path 'delta_scatter.jpeg' , emit: delta_scatter , optional: true

     """
       delta_plot.py  -t delta_PB_two.csv  -b delta_baci.csv  -c delta_CFML.csv  -g delta_Gubbins.csv
     """
}


process TreeCmp_summary {

     publishDir "${PWD}/Summary_${params.outDir}", mode: "copy"
     maxForks 1


     input:
        path PB_dist
        path Dist_Gubbins
        path Dist_CFML
        path PB_wrf
        path Gubbins_wrf
        path CFML_wrf



     output:
        path   'Dist_summary.jpeg'    , emit: FigTreeDist

     """
       cmpTree_plot.py  -p dist_PB.csv  -g dist_Gubbins.csv  -m dist_CFML.csv  -prf wrf_PB.csv  -grf  wrf_Gubbins.csv -mrf wrf_CFML.csv
     """
}

workflow Sim {
        take:
            iteration
            frequencies
            rates
            nu_sim
            recomlen
            recomrate
        main:
            MakeClonalTree(iteration)
            BaciSim(MakeClonalTree.out.Clonaltree,nu_sim,recomlen,recomrate)
            Seq_gen(BaciSim.out.BaciSimtrees,frequencies,rates,BaciSim.out.nu_sim,BaciSim.out.recomlen,BaciSim.out.recomrate)

        emit:
            clonaltree = BaciSim.out.Clonaltree
            genome = Seq_gen.out.Wholegenome
            recom_log = BaciSim.out.Recomlog
            iteration = BaciSim.out.Range
            recom_stat = BaciSim.out.Recomstat
            nu_sim = BaciSim.out.nu_sim
            recomlen = BaciSim.out.recomlen
            recomrate= BaciSim.out.recomrate
            unroot_clonaltree = BaciSim.out.unroot_Clonaltree

}


workflow Gubbins {
        take:
            genome
            clonaltree
            recom_log
            iteration
            nu_sim
            recomlen
            recomrate
        main:
            Run_Gubbins(genome,iteration,nu_sim,recomlen,recomrate)
            Gubbins_result(clonaltree,recom_log,genome,Run_Gubbins.out.Gubbinstree,Run_Gubbins.out.GubbinsRecom,Run_Gubbins.out.GubbinsStat,nu_sim,recomlen,recomrate)
        emit:
            RMSE_Gubbins = Gubbins_result.out.Rmse_Gubbins
            GubbinsRescaletree = Gubbins_result.out.GubbinsRescaletree
            Rcount_Gubbins = Gubbins_result.out.Gubb_rcount
            Delta_Gubbins = Gubbins_result.out.Gubb_delta
            Dist_Gubbins = Gubbins_result.out.Gubb_dist
            WRF_Gubbins = Gubbins_result.out.Gubb_wrf

}


workflow ClonalFrameML {
        take:
            clonaltree
            recom_log
            genome
            raxml_tree
            iteration
            nu_sim
            recomlen
            recomrate

        main:
            CFML(genome,raxml_tree,iteration,nu_sim,recomlen,recomrate)
            CFML_result(genome,CFML.out.CFML_recom,CFML.out.CFMLtree,recom_log,clonaltree,iteration,nu_sim,recomlen,recomrate)
        emit:
            CFMLtree = CFML.out.CFMLtree
            RMSE_CFML = CFML_result.out.RMSE_CFML
            Rcount_CFML = CFML_result.out.CFML_rcount
            Delta_CFML = CFML_result.out.CFML_delta
            Dist_CFML = CFML_result.out.CFML_dist
            WRF_CFML = CFML_result.out.CFML_wrf
}



workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }
    if (params.mode == 'sim') {
        println "Detect recombination in simulated data..."
        params.simulation = 1
        Sim(iteration,frequencies,rates,nu_sim,recomlen,recomrate)
        Get_raxml_tree(Sim.out.genome,Sim.out.iteration,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)

        if (params.method =~ /cfml/) {
            ClonalFrameML(Sim.out.clonaltree,Sim.out.recom_log,Sim.out.genome,Get_raxml_tree.out.MyRaxML,Sim.out.iteration,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)
            CollectedRMSE_CFML = ClonalFrameML.out.RMSE_CFML.collectFile(name:"rmse_CFML.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedRcount_CFML = ClonalFrameML.out.Rcount_CFML.collectFile(name:"rcount_CFML.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedDelta_CFML = ClonalFrameML.out.Delta_CFML.collectFile(name:"delta_CFML.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:true , sort: false)
            CollectedDist_CFML = ClonalFrameML.out.Dist_CFML.collectFile(name:"dist_CFML.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedWRF_CFML = ClonalFrameML.out.WRF_CFML.collectFile(name:"wrf_CFML.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)

        }
        if (params.method =~ /gub/) {
            Gubbins(Sim.out.genome,Sim.out.clonaltree,Sim.out.recom_log,Sim.out.iteration,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate)
            CollectedRMSE_Gubb = Gubbins.out.RMSE_Gubbins.collectFile(name:"rmse_Gubbins.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedRcount_Gubb = Gubbins.out.Rcount_Gubbins.collectFile(name:"rcount_Gubbins.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedDelta_Gubb = Gubbins.out.Delta_Gubbins.collectFile(name:"delta_Gubbins.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:true , sort: false)
            CollectedDist_Gubb = Gubbins.out.Dist_Gubbins.collectFile(name:"dist_Gubbins.csv",storeDir:"${PWD}/Summary_${params.outDir}",keepHeader:false,sort:false)
            CollectedWRF_Gubb = Gubbins.out.WRF_Gubbins.collectFile(name:"wrf_Gubbins.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
        }
        if (params.method =~ /pb/) {
            PhiloBacteria(Sim.out.clonaltree,Sim.out.recom_log,Sim.out.genome,Get_raxml_tree.out.MyRaxML,Sim.out.nu_sim,Sim.out.recomlen,Sim.out.recomrate,Sim.out.recom_stat)
            CollectedRMSE_PB_two = PhiloBacteria.out.PB_RMSE_two.collectFile(name:"rmse_PB_two.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedRcount_Baci = PhiloBacteria.out.Baci_rcount.collectFile(name:"rcount_baci.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedRcount_PB_two = PhiloBacteria.out.PB_rcount_two.collectFile(name:"rcount_PB_two.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedDelta_Baci = PhiloBacteria.out.Baci_Delta.collectFile(name:"delta_baci.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:true , sort: false)
            CollectedDelta_PB_two = PhiloBacteria.out.PB_Delta_two.collectFile(name:"delta_PB_two.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:true , sort: false)
            CollectedDist_PB = PhiloBacteria.out.PB_dist.collectFile(name:"dist_PB.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)
            CollectedWRF_PB = PhiloBacteria.out.PB_wrf.collectFile(name:"wrf_PB.csv",storeDir:"${PWD}/Summary_${params.outDir}", keepHeader:false , sort: false)

        }

         if (params.analyse == true)  {
            RMSE_summary(CollectedRMSE_CFML,CollectedRMSE_Gubb,CollectedRMSE_PB_two)
            Delta_Summary(CollectedDelta_Baci,CollectedDelta_PB_two,CollectedDelta_Gubb,CollectedDelta_CFML)
            TreeCmp_summary(CollectedDist_PB,CollectedDist_Gubb,CollectedDist_CFML,CollectedWRF_PB,CollectedWRF_Gubb,CollectedWRF_CFML)
         }
    }
     if (params.mode == 'emp') {
        println "Detect recombination in empirical sequence alignments..."
        params.simulation = 0
        recom_log = tuple(1,dummy_file)
        clonaltree = dummy_file
        genome = params.seq
        recomstat = tuple(1,dummy_file)

        Get_raxml_tree(genome,1,1,1,1)

        if (params.method =~ /cfml/) {
            //ClonalFrameML(clonaltree,recom_log,genome,Get_raxml_tree.out.MyRaxML,1,1,1,1)
            CFML(genome,Get_raxml_tree.out.MyRaxML,1,1,1,1)
        }
        if (params.method =~ /gub/) {
            Gubbins(genome,clonaltree,recom_log,1,1,1,1)
        }
        if (params.method =~ /pb/) {
            PhiloBacteria(clonaltree,recom_log,genome,Get_raxml_tree.out.MyRaxML,1,1,1,recomstat)
        }

    }

}