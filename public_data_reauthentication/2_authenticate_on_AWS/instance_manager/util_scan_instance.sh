aws ec2 describe-instances \
  --filters "Name=instance-state-name,Values=running" \
  --query 'Reservations[*].Instances[*].[InstanceId, PublicIpAddress, InstanceType, LaunchTime, Tags[?Key==`Name`].Value | [0]]' \
  --output text 2>/dev/null | sort -k4,4V
