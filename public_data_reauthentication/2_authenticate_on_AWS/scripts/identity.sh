TOKEN=`curl -s -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600"`
INSTANCEID=`curl -s -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id`
INSTANCENAME=`aws ec2 describe-instances --instance-ids $INSTANCEID --query 'Reservations[0].Instances[0].Tags[?Key==\`Name\`].Value' --output text`
BATCH=`echo $INSTANCENAME | cut -d"_" -f2`
INDEX=`echo $INSTANCENAME | cut -d"_" -f3`
