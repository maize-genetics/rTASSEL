# Build from rocker's tidyverse stack (v4.5)
FROM rocker/tidyverse:4.5

# rTASSEL Git tag (e.g. v0.12.0); passed by CI at image build time
ARG RTASSEL_TAG=v0.12.0

# Install system libraries for Java, rJava, Git, and XML/SSL support
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      default-jdk \
      git \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev && \
    rm -rf /var/lib/apt/lists/*

# Ensure JAVA_HOME is set for rJava; align R's Java config with CI (setup-tassel action)
ENV JAVA_HOME=/usr/lib/jvm/default-java
RUN R CMD javareconf

# Install pak
RUN R -e "install.packages('pak')"

# Install rTASSEL (RTASSEL_TAG expanded by the shell, e.g. v0.12.0)
RUN R -e "pak::pak('maize-genetics/rTASSEL@${RTASSEL_TAG}')"

# Setup rTASSEL
RUN R -e "rTASSEL::setupTASSEL()"
