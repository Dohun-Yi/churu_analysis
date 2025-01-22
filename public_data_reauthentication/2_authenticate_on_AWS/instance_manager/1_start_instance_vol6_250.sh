#NOTE: n_instance must be synced with SRA file split (0_prep~~sh) and the printf below (index=~~)!! bad code, but won't fix it
n_instance=40
batch=$(date +%s)

# Using throoughput 250MiB/sec and 3000 IOPS.
launch_instance(){ # f g h i j k
    volumn1="{ \"DeviceName\": \"/dev/sdf\", \"Ebs\": { \"VolumeSize\": 200, \"VolumeType\": \"gp3\", \"Throughput\": 250, \"Iops\": 3000 } }" # gp3, .2 TB 250MiB/s, 3000 IOPS
    volumn2="{ \"DeviceName\": \"/dev/sdg\", \"Ebs\": { \"VolumeSize\": 200, \"VolumeType\": \"gp3\", \"Throughput\": 250, \"Iops\": 3000 } }" # gp3, .2 TB 250MiB/s, 3000 IOPS
    volumn3="{ \"DeviceName\": \"/dev/sdh\", \"Ebs\": { \"VolumeSize\": 200, \"VolumeType\": \"gp3\", \"Throughput\": 250, \"Iops\": 3000 } }" # gp3, .2 TB 250MiB/s, 3000 IOPS
    volumn4="{ \"DeviceName\": \"/dev/sdi\", \"Ebs\": { \"VolumeSize\": 200, \"VolumeType\": \"gp3\", \"Throughput\": 250, \"Iops\": 3000 } }" # gp3, .2 TB 250MiB/s, 3000 IOPS
    volumn5="{ \"DeviceName\": \"/dev/sdj\", \"Ebs\": { \"VolumeSize\": 200, \"VolumeType\": \"gp3\", \"Throughput\": 250, \"Iops\": 3000 } }" # gp3, .2 TB 250MiB/s, 3000 IOPS
    volumn6="{ \"DeviceName\": \"/dev/sdk\", \"Ebs\": { \"VolumeSize\": 200, \"VolumeType\": \"gp3\", \"Throughput\": 250, \"Iops\": 3000 } }" # gp3, .2 TB 250MiB/s, 3000 IOPS
    volumns="[$volumn1, $volumn2, $volumn3, $volumn4, $volumn5, $volumn6]"
    instance_type="c6g.8xlarge"

    batch=$1
    index=$(printf "%02d" $2) ### here
    instance_name="ChuruPub_${batch}_${index}"

    echo "Launching $instance_name..."
    aws ec2 run-instances \
        --image-id ami-0fdcbfc2802f642d3 \
        --instance-type $instance_type \
        --key-name etching \
        --security-groups churu-group \
        --iam-instance-profile Name=churu-reader \
        --block-device-mappings "$volumns" \
        --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$instance_name}]" \
        --user-data \
"#!/bin/bash
export NAME=$instance_name
export NUMBER=$index

POST_DATA()
{
  cat <<EOF
{
  \"text\":\"\$1\"
}
EOF
}
notify() {
    MESSAGE=\"\${NAME}_\$1\"
    # Send message to my Slack
    # curl -X POST -H 'Content-type: application/json' --data \"\$(POST_DATA \$MESSAGE)\" {URL}
}
notify hello

## mount
mkfs -t xfs /dev/sdf &
mkfs -t xfs /dev/sdg &
mkfs -t xfs /dev/sdh &
mkfs -t xfs /dev/sdi &
mkfs -t xfs /dev/sdj &
mkfs -t xfs /dev/sdk &
wait
mkdir /data1
mkdir /data2
mkdir /data3
mkdir /data4
mkdir /data5
mkdir /data6
mount /dev/sdf /data1
mount /dev/sdg /data2
mount /dev/sdh /data3
mount /dev/sdi /data4
mount /dev/sdj /data5
mount /dev/sdk /data6
cd /data1
notify format_done


# install
aws s3 cp --recursive s3://dohun/churu/setup /data1/mock/setup
aws s3 cp --recursive s3://dohun/instance_resource /data1/instance_resource/
aws s3 cp --recursive s3://dohun/churu /data1/churu &
cd ./mock/setup/
ln -s /data1/instance_resource/scripts/setup_sratools.sh
bash ./install_yum.sh
bash ./install_pip.sh &
bash ./install_bwa.sh &
bash ./install_star.sh &
bash ./install_kmc.sh &
bash ./install_samtools.sh &
bash ./setup_sratools.sh
wait
ln -srf /data1/mock/* /data1/churu/
chmod +x /data1/churu/churu.sh
export PATH=\$PATH:/usr/local/ncbi/sra-tools/bin/
cd -
notify setup_done


ln -srf ./instance_resource/SRA_Accessions.txt_split$index ./instance_resource/SRA_Accessions.batch.txt

mkdir -p ./log
bash ./instance_resource/scripts/master.sh 1> ./log/master_back.out 2> ./log/master_back.err &
job1=\$!
sleep 5
bash ./instance_resource/scripts/logger.sh 1> ./log/log.out 2> ./log/log.err &
job2=\$!
notify started_all_processor

sar 5 10000000 1> ./log/sar.out &
wait \$job1
kill $job2
bash ./instance_resource/scripts/logger_once.sh
notify i_think_im_done
bash ./instance_resource/scripts/suicide.sh
"
}

export -f launch_instance
seq 1 $n_instance | parallel -j $n_instance launch_instance ${batch}
