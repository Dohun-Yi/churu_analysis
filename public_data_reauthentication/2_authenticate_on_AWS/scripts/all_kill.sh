pgrep -f 'downloader.sh' | xargs kill
pgrep -f 'decompressor.sh' | xargs kill
pgrep -f 'runner.sh' | xargs kill
pgrep -f 'logger.sh' | xargs kill
pgrep -f 'fasterq-dump' | xargs kill -9
pgrep -f 'STAR' | xargs kill -9
rm -rf /data/downloading/
rm -rf /data/downloaded/
rm -rf /data/decomped/
rm -rf /data/processing/
rm -rf /data/log/
rm -rf /data/tmp/
