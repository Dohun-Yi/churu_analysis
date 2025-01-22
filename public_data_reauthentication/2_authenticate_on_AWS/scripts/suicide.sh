source $(dirname $0)/identity.sh
aws ec2 terminate-instances --instance-ids $INSTANCEID
