#!/bin/sh
sudo mkdir -p /mnt/cinder/$1;
sudo mkfs.xfs /dev/vdb;
sudo mount /dev/vdb /mnt/cinder/$1;
sudo mkswap /mnt/cinder/$1/swapfile;
sudo swapon /mnt/cinder/$1/swapfile;
