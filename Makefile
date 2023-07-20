init:
	git lfs install
	git lfs pull
	conda install git pip
	pip install -r requirements.txt
	if [ ! -d TP_prediction/output ]; then mkdir TP_prediction/output; fi