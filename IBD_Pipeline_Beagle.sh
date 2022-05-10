#This is a pipeline I modified and ran to call Identical-By-Descent tracts in eastern Africa SNP array data using the program Beagle
#The original pipeline and associated scripts can be found here: https://github.com/halasadi/ibd_data_pipeline/blob/master/runibdpipeline
#Ultimately I used hap-ibd since it produced calls that were more accurate for my dataset

# chromosome 1-22
CHROMOSOMES = [str(i) for i in range(1,23)]

# testing
#CHROMOSOMES = [str(i) for i in range(20,21)]

SEEDS = ['1', '5', '10', '15', '20']

### CONFIG FILE ###
configfile: "config.json"
# 'ibdcalls' specifies the directory where the ibdcalls are located
# 'prefix' specifies the prefix the script will append to the output
# 'plinkfile' specifies the location of plink files (broken up by chr) and the prefix of the plink files
# 'build' the genome build of the data (b36 or b37)

rule all:
    input:
        expand(config['ibdcalls'] + '/' + config['prefix'] + '_{chr}.merged.nosparse_cm_finalqc.ibd', chr=CHROMOSOMES)

rule remove_gaps:
    input:
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}.merged.nosparse_cm.ibd'
    output:
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}.merged.nosparse_cm_finalqc.ibd'
    params:
        gapthres = '0.4'
    shell:
        'python remove-gaps-fibd.py -b {input[0]} -o {output[0]} -g {params.gapthres}'

rule convert_to_cM:
    input:
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}.merged.nosparse.ibd'
    output:
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}.merged.nosparse_cm.ibd'
    params:
        chromo='{chr}'
    shell:
        'Rscript interpolate.R {params.chromo} {config[build]} {input} {output}'

rule remove_sparse:
    input:
        config['ibdcalls'] + '/' + 'chr{chr}.sparseregions', config['ibdcalls'] + '/' + config['prefix'] + '_{chr}.merged.ibd'
    output:
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}.merged.nosparse.ibd'
    shell:
        'python remove_sparse.py -s {input[0]} -i {input[1]} > {output[0]}'

rule create_sparseregions:
    input:
        config['plinkDir'] + '/broken_up_by_chr/' + config['plinkPrefix'] +  '_{chr}.bim'
    output:
        config['ibdcalls'] + '/' + 'chr{chr}.sparseregions'
    shell:
        'python sparse_regions.py -b {input[0]} > {output[0]}'

rule merge_ibd_segments:
    input:
        config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[0] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[0] + '.ibd',
        config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[1] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[1] + '.ibd',
        config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[2] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[2] + '.ibd'
        config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[3] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[3] + '.ibd'
        config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[4] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[4] + '.ibd'

    output:
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}.merged.ibd'
    params:
        ibdmerge = '/share/hennlab/progs/ibdmerge.jar' #original script calls for ibdmerge.jar
    shell:
        'cat {input} | java -jar {params.ibdmerge} > {output}'

rule relabelIBD:
    input:
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[0] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[0] + '.ibd',
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[1] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[1] + '.ibd',
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[2] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[2] + '.ibd'
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[3] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[3] + '.ibd'
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[4] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[4] + '.ibd'
    output:
        config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[0] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[0] + '.ibd',
        config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[1] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[1] + '.ibd',
        config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[2] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[2] + '.ibd'
        config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[3] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[3] + '.ibd'
        config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[4] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '.relabeled.{chr}_' + SEEDS[4] + '.ibd'
    shell:
        'python rename.py {input[0]} > {output[0]} && '
        'python rename.py {input[1]} > {output[1]} && '
        'python rename.py {input[2]} > {output[2]} && '
        'python rename.py {input[3]} > {output[3]} && '
        'python rename.py {input[4]} > {output[4]} && '
        'python rename.py {input[5]} > {output[5]} && '
        'python rename.py {input[6]} > {output[6]} && '
        'python rename.py {input[7]} > {output[7]} && '
        'python rename.py {input[8]} > {output[8]} && '
        'python rename.py {input[9]} > {output[9]} '

rule callIBD:
    input:
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}.vcf'
    output:
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[0] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[0] + '.ibd',
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[1] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[1] + '.ibd',
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[2] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[2] + '.ibd'
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[3] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[3] + '.ibd'
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[4] + '.hbd', config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_' + SEEDS[4] + '.ibd'
    params:
        f = config['ibdcalls'] + '/' + config['prefix'] + '_{chr}_',
        heap = '22',
        beagle = '/share/hennlab/progs/beagle.27Jan18.7e1.jar',
	nthreads = '4',
    shell:
        'java -Xmx{params.heap}g -jar {params.beagle} gt={input} ibd=true seed={SEEDS[0]} nthreads={params.nthreads} out={params.f}{SEEDS[0]} && '
        'java -Xmx{params.heap}g -jar {params.beagle} gt={input} ibd=true seed={SEEDS[1]} nthreads={params.nthreads} out={params.f}{SEEDS[1]} && '
        'java -Xmx{params.heap}g -jar {params.beagle} gt={input} ibd=true seed={SEEDS[2]} nthreads={params.nthreads} out={params.f}{SEEDS[2]} && '
        'java -Xmx{params.heap}g -jar {params.beagle} gt={input} ibd=true seed={SEEDS[3]} nthreads={params.nthreads} out={params.f}{SEEDS[3]} && '
        'java -Xmx{params.heap}g -jar {params.beagle} gt={input} ibd=true seed={SEEDS[4]} nthreads={params.nthreads} out={params.f}{SEEDS[4]} && '
        'gunzip {output[0]}.gz && '
        'gunzip {output[1]}.gz && '
        'gunzip {output[2]}.gz && '
        'gunzip {output[3]}.gz && '
        'gunzip {output[4]}.gz && '
        'gunzip {output[5]}.gz && '
        'gunzip {output[6]}.gz && '
        'gunzip {output[7]}.gz && '
        'gunzip {output[8]}.gz && '
        'gunzip {output[9]}.gz '

rule convert_to_vcf:
    input:
        config['plinkDir'] + '/' + config['plinkPrefix'] + '.bed',
        config['plinkDir'] + '/' + config['plinkPrefix'] + '.bim',
        config['plinkDir'] + '/' + config['plinkPrefix'] + '.fam'
    output:
        config['ibdcalls'] + '/' + config['prefix'] + '_{chr}.vcf'
    params:
        chromo = '{chr}',
        plink = 'plink',
        out = config['ibdcalls'] + '/' + config['prefix'] + '_{chr}'
    shell:
        '{params.plink} --bfile ' + config['plinkDir'] + '/' + config['plinkPrefix'] + '  --chr {params.chromo} --recode vcf --out {params.out}'

rule break_up_plink:
    input:
        config['plinkDir'] + '/' + config['plinkPrefix'] + '.bed',
        config['plinkDir'] + '/' + config['plinkPrefix'] + '.bim',
        config['plinkDir'] + '/' + config['plinkPrefix'] + '.fam'
    output:
        config['plinkDir'] + '/broken_up_by_chr/' + config['plinkPrefix'] +  '_{chr}.bim'
    params:
        chromo = '{chr}',
        plink = 'plink',
        out = config['plinkDir'] + '/broken_up_by_chr/' + config['plinkPrefix'] +  '_{chr}'
    shell:
        '{params.plink} --bfile ' + config['plinkDir'] + '/' + config['plinkPrefix'] + ' --chr {params.chromo} --make-bed --out {params.out}'
