#!/bin/bash

# 独立提交脚本：无需依赖 a.sh / b.sh / c.sh。
# 作用：
# 1) 自动创建 a、b、c 三个目录
# 2) 在各自目录内生成与原脚本一致的 L*/h* 目录结构、paramC_sets.txt、dqmc，并提交任务
# 3) 默认依次执行 a、b、c；也可传参只执行某些集合。例如：
#    ./set.sh a c

# set script_dir (robust even when sourced)
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

set -e

# 全局通用参数（与原脚本保持一致）
# 仅跑 L=12 尺寸
# L_list=(14)
L_list=(10)
# beta_list=(32)
beta_list=(28)
NWRAP=5
nbin=7500
nsweep=1
Nthermal=0
Nwarm=500
Nglobal=8
absolute=5
sector_charge=1
max_cpus_per_node=4  # 每节点CPU上限

# 三组 h 列表（对应 a/b/c）
# Nat Phys 2017 Fig. 3: 扫 h 范围 [0.52, 0.62]，临界场 h_c≈0.56(2)
# declare -a h_list_a=($(seq 0.52 0.02 0.62))
declare -a h_list_a=( 0.52 )
# declare -a h_list_a=($(seq 2.5 0.3 4.0))
declare -a h_list_b=(1.5 1.75 2)
declare -a h_list_c=(1.79985 1.89995 2 2.1016 2.2018 2.30055 2.40075 2.50095 2.6012)

# 队列选择相关函数（基于 pbsnodes 与 qstat -f，兼容你们的 pbsfree.sh 统计口径）
# 支持通过环境变量覆盖节点属性与队列名：TESTQ_NODE_PROPERTY/WORKQ_NODE_PROPERTY，TESTQ_NAME/WORKQ_NAME

# 使用 pbsnodes 聚合计算每个队列的空闲CPU核心数
get_queue_free_cpus() {
    local queue="$1"
    local TESTQ_NODE_PROPERTY="${TESTQ_NODE_PROPERTY:-testq}"
    local WORKQ_NODE_PROPERTY="${WORKQ_NODE_PROPERTY:-workq}"
    local TESTQ_NAME="${TESTQ_NAME:-testq}"
    local WORKQ_NAME="${WORKQ_NAME:-workq}"

    if ! command -v pbsnodes >/dev/null 2>&1; then
        echo ""; return 1
    fi

    pbsnodes -a 2>/dev/null | awk -v target="$queue" -v test_prop="$TESTQ_NODE_PROPERTY" -v work_prop="$WORKQ_NODE_PROPERTY" -v test_name="$TESTQ_NAME" -v work_name="$WORKQ_NAME" '
    BEGIN {}
    /^[a-zA-Z0-9]/ { node=$1; nodes[node]=1 }
    /state =/ { node_state[node]=$3 }
    /resources_available.ncpus =/ { node_avail[node]=$3+0 }
    /resources_assigned.ncpus =/ { node_assigned[node]=$3+0 }
    /properties =/ {
        props_line=$0; sub(/^.*properties =[ ]*/, "", props_line); node_props[node]=props_line
    }
    /queue =/ { node_queue[node]=$3 }
    END {
        for (n in nodes) {
            avail=(n in node_avail?node_avail[n]:0)
            assigned=(n in node_assigned?node_assigned[n]:0)
            props=(n in node_props?node_props[n]:"")
            props_cmp=props; gsub(/[[:space:]]+/, "", props_cmp)
            q_label=""
            if (n in node_queue) {
                qv=node_queue[n]
                if (tolower(qv)==tolower(test_name)) q_label="testq";
                else if (tolower(qv)==tolower(work_name)) q_label="workq";
            }
            if (q_label=="") {
                if (test_prop!="" && props_cmp ~ ("(^|,)" test_prop "(,|$)")) q_label="testq";
                else if (work_prop!="" && props_cmp ~ ("(^|,)" work_prop "(,|$)")) q_label="workq";
            }
            if (q_label=="") q_label="workq";
            queue_total[q_label]+=avail;
            queue_used[q_label]+=assigned;
        }
        total_work=("workq" in queue_total?queue_total["workq"]:0)
        used_work=("workq" in queue_used?queue_used["workq"]:0)
        total_test=("testq" in queue_total?queue_total["testq"]:0)
        used_test=("testq" in queue_used?queue_used["testq"]:0)
        free_work=total_work-used_work; if (free_work<0) free_work=0
        free_test=total_test-used_test; if (free_test<0) free_test=0
        if (tolower(target)=="workq") print free_work; else if (tolower(target)=="testq") print free_test; else print 0
    }'
}

# 使用 qstat -f 统计队列中的排队作业数量（Q 或 H）
get_queue_queued_jobs() {
    local queue="$1"
    if ! command -v qstat >/dev/null 2>&1; then
        echo 0; return 0
    fi
    qstat -f 2>/dev/null | awk -v target="$queue" '
    BEGIN{q=""; s=""; c=0}
    /^Job Id:/ { q=""; s="" }
    /queue =/ { q=$3 }
    /job_state =/ { s=$3; if (q==target && (s=="Q" || s=="H")) c++ }
    END{ print c+0 }'
}

# 选择最佳队列：
# 1) 若仅一方存在空闲CPU(>0)，选该队列
# 2) 若两者都无空闲或都存在空闲，则比较排队更少者
# 3) 若仍相同，按空闲CPU更多者；再相同默认选 workq
choose_queue() {
    local q1="workq" q2="testq"
    local free1 free2 queued1 queued2
    free1=$(get_queue_free_cpus "$q1"); free2=$(get_queue_free_cpus "$q2")
    # 非数字时置0
    [[ "$free1" =~ ^[0-9]+$ ]] || free1=0
    [[ "$free2" =~ ^[0-9]+$ ]] || free2=0

    if (( free1>0 && free2==0 )); then echo "$q1"; return 0; fi
    if (( free2>0 && free1==0 )); then echo "$q2"; return 0; fi

    queued1=$(get_queue_queued_jobs "$q1"); queued2=$(get_queue_queued_jobs "$q2")
    [[ "$queued1" =~ ^[0-9]+$ ]] || queued1=0
    [[ "$queued2" =~ ^[0-9]+$ ]] || queued2=0

    if (( queued1 < queued2 )); then echo "$q1"; return 0; fi
    if (( queued2 < queued1 )); then echo "$q2"; return 0; fi

    if (( free1 > free2 )); then echo "$q1"; return 0; fi
    if (( free2 > free1 )); then echo "$q2"; return 0; fi

    echo "$q1"
}

# 预检：选择可执行文件，并检查依赖文件是否存在
bin_candidates=(z2_gauge.out dqmc_gauge)
bin_name=""
for cand in "${bin_candidates[@]}"; do
    if [[ -f "${script_dir}/${cand}" ]]; then
        bin_name="$cand"
        break
    fi
done
if [[ -z "${bin_name}" ]]; then
    echo "错误: 未在 ${script_dir} 找到可执行文件 dqmc_gauge。请将编译好的可执行文件放到 ${script_dir} 后重试。" >&2
    exit 1
fi
for f in seeds.txt confin.txt; do
    if [[ ! -f "${script_dir}/${f}" ]]; then
        echo "错误: 未找到文件 ${script_dir}/${f}。请将该文件放到脚本所在目录后重试。" >&2
        exit 1
    fi
done

# 解析要运行的集合：默认仅运行 a 组
targets=("a")
if [[ $# -gt 0 ]]; then
    targets=("$@")
fi

# 生成与提交函数
run_group() {
    local tag="$1"            # 1 / 2 / 3
    local job_suffix="$2"     # 1 / 2 / 3（用于队列名）
    local Jg_value="$3"       # 第一行中的 Jg
    local mu_value="$4"       # 第一行中的 mu
    shift 4
    local group_h_list=("$@")  # 剩余参数为该组 h 列表

    local subdir="${script_dir}/${tag}"
    mkdir -p "${subdir}"
    (
        cd "${subdir}"

        RT="${RT:-0.5}"
        for ((i=0; i<${#L_list[@]}; i++)); do
            L=${L_list[$i]}
            cpu=1

            mkdir -p L${L}
            cd L${L}

            for ((jj=0; jj<${#group_h_list[@]}; jj++)); do
                printf -v h_fmt "%.2f" "${group_h_list[$jj]}"
                mkdir -p h${h_fmt}
                cd h${h_fmt}

                # pre calc
                h="${h_fmt}"

                beta=${beta_list[$i]}
                Ltrot=$(awk -v b="$beta" 'BEGIN{ printf "%.0f", 20*b }')

                # 计算所需的节点数和每个节点的任务数
                nodes=$(( (cpu + max_cpus_per_node - 1) / max_cpus_per_node ))
                ntasks_per_node=$(( (cpu + nodes - 1) / nodes ))

                # paramC_sets（注意将 ThetaX ThetaY 行放在注释之前）
                cat > paramC_sets.txt << EOF
${RT}  ${Jg_value} ${h}        ${mu_value}   ${sector_charge}
${L}   ${L}        ${Ltrot}    ${beta}
${L}   ${L}        ${Ltrot}
${NWRAP}           ${nbin}     ${nsweep}
.false.            ${Nthermal}
.false.            ${Nwarm}
.true.             ${Nglobal}
${absolute}

# RT, J, h, mu, Q
# Nlx, Nly, Ltrot, Beta
# NlxTherm, NlyTherm, LtrotTherm
# Nwrap, Nbin, Nsweep
# is_tau, Nthermal
# is_warm, Nwarm
# is_global, Nglobal
# absolute
EOF

                # dqmc
                temp_NPROCS='`wc -l < \$PBS_NODEFILE` '
                selected_queue="$(choose_queue)"
                # selected_queue="testq"
                # selected_queue="workq"
                cat > dqmc << EOF
#!/bin/bash
#PBS -N z2gauge2-1_L${L}h${h}beta${beta}
#PBS -l nodes=${nodes}:ppn=${ntasks_per_node}
#PBS -q ${selected_queue}
#PBS -l walltime=2400:00:00
#PBS -j oe
cd \$PBS_O_WORKDIR
ulimit -s unlimited
ulimit -c unlimited
NPROCS=${temp_NPROCS}
mpirun -machinefile \$PBS_NODEFILE -np \$NPROCS ./${bin_name}
EOF

                # copy 依赖与可执行文件
                cp -- "${script_dir}/${bin_name}" "${script_dir}/seeds.txt" "${script_dir}/confin.txt" ./

                # 提交
                pwd
                qsub dqmc

                cd ../
            done

            cd ../
        done
    )
}

for t in "${targets[@]}"; do
    case "${t}" in
        a)
            run_group "a" "a" "1.0" "0.3" "${h_list_a[@]}" ;;
        b)
            # 若用户显式传入 b/c，则仍支持原有运行
            run_group "b" "b" "0.5" "1.5" "-2.0" "${h_list_b[@]}" ;;
        c)
            run_group "c" "c" "0.5" "3.0" "-2.0" "${h_list_c[@]}" ;;
        *)
            echo "警告: 忽略未知目标 '${t}'（仅支持 a/b/c）" >&2 ;;
    esac
done

echo "全部完成。"
