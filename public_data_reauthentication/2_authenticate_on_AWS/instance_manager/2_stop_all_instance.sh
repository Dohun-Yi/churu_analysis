#!/bin/bash

# Get a list of all running instances
instance_ids=$(aws ec2 describe-instances \
  --query 'Reservations[].Instances[?State.Name==`running`].InstanceId' \
  --output text)

# Loop through the instance IDs and stop each instance
for instance_id in $instance_ids; do
  echo "terminating $instance_id"
  aws ec2 terminate-instances --instance-ids "$instance_id"
  # aws ec2 stop-instances --instance-ids "$instance_id"
done

