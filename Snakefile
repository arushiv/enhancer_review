import pandas
import os

DIRECTORIES = {
    'intermediateFiles': "intermediateFiles",
    'figures' : "figures",
    }

SCRIPTS = {
    'plot' : "scripts/heatmap.R"
}

CELLS = ['GM12878', 'H1', 'HepG2', 'HMEC', 'HSMM', 'Huvec', 'K562', 'NHLF']
ANNOTATIONS = ['Stretch_Enhancers', 'Super_Enhancers', 'Typical_Enhancers']

rule final:
    input:
        os.path.join(DIRECTORIES['figures'], "fig.pdf")

rule getFractions:
    input:
        annot1 = "/lab/work/arushiv/annotations/{cell1}.{annotation1}.bed",
        annot2 = "/lab/work/arushiv/annotations/{cell2}.{annotation2}.bed"
    output:
        main = temp(os.path.join(DIRECTORIES['intermediateFiles'], "{cell1}.{annotation1}.intersect.{cell2}.{annotation2}.bed"))
    shell:
        r"""
        coverage=`intersectBed -a {annot1} -b {annot2} -wao | awk '{tot=tot+$7}END{print tot}'`;
        length=`less {annot1} | awk '{len=len+$3-$2}END{print len}'` ;
        region_overlap=`intersectBed -a {annot1} -b {annot2} -u | wc -l` ;
        total_regions=`less {annot1} | wc -l` ;
        echo -e "{cell1}\t{annotation1}\t{cell2}\t{annotation2}\t$length\t$coverage\t${{total_regions}}\t${{region_overlap}}" > {output.main}
        """

rule mergeFractions:
    input:
        expand(rules.getFractions.output.main,
               cell1 = CELLS,
               cell2 = CELLS,
               annotation1 = ANNOTATIONS,
               annotation2 = ANNOTATIONS)
    output:
        tempfile = temp(os.path.join(DIRECTORIES['intermediateFiles'], "coverageFraction_all.dat.temp")),
        main = os.path.join(DIRECTORIES['intermediateFiles'], "coverageFraction_all.dat")
    shell:
        r"""
        cat {input} | awk '{{if (( $1==$3 && $2==$4 )) {{$6=$8="NA"}}; print $0}}' OFS='\t' > {output.tempfile};
        echo -e "Cell_type1\tEnhancerType1\tCell_type2\tEnhancerType2\ttotal_length1\tbase_overlap\ttotal_regions1\tregion_overlap" | cat - {output.tempfile} > {output.main}
        """

rule plot:
    input:
        rules.mergeFractions.output.main,
    output:
        os.path.join(DIRECTORIES['figures'], "fig.pdf")
    params:
        script = SCRIPTS['plot']
    shell:
        """
        Rscript {params.script} {input} {output}
        """
        
