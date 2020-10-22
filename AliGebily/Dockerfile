# Dockerized version of DECONV > DE transformation code, provided by Harvard.
FROM shotat/anaconda3-xgboost-docker

RUN pip install --upgrade xgboost
COPY ./* /python-model/
ADD ./resources/ /python-model/resources/
# Make the file readable, writable and executable by everyone
RUN chmod 777 /python-model/predict.sh

ENTRYPOINT [ "/python-model/predict.sh" ]
