pipeline {
  agent {
    dockerfile {
      additionalBuildArgs '--tag akantu-environment'
    }
    
  }
  stages {
    stage('Configure') {
      steps {
        sh 'env'
        sh 'mkdir -p build'
        sh 'cd build; cmake -DAKANTU_COHESIVE_ELEMENT:BOOL=TRUE -DAKANTU_IMPLICIT:BOOL=TRUE -DAKANTU_PARALLEL:BOOL=TRUE -DAKANTU_PYTHON_INTERFACE:BOOL=TRUE -DAKANTU_TESTS:BOOL=TRUE ..'
      }
    }
    stage('Compile') {
      steps {
	wrap([$class: 'AnsiColorBuildWrapper', 'colorMapName': 'xterm']) {
          sh 'make -C build/src || true'
	}
      }
    }

    stage ('Warnings gcc') {
      steps {
        warnings(consoleParserss: [[parserName: 'GNU Make + GNU C Compiler (gcc)']])
      }
    }

    stage('Compile python') {
      steps {
        sh 'make -C build/python || true'
      }
    }

    stage('Compile tests') {
      steps {
        sh 'make -C build/test || true'
      }
    }

    stage('Tests') {
      steps {
          sh 'curl https://raw.githubusercontent.com/rpavlik/jenkins-ctest-plugin/master/ctest-to-junit.xsl -o ctest-to-junit.xsl'
	  sh 'cd build/ && ctest -T test --no-compress-output || true'
	  sh 'xsltproc ctest-to-junit.xsl  build/Testing/`head -n 1 < build/Testing/TAG`/Test.xml > CTestResults.xml'
      }
    }
  }
  environment {
    BLA_VENDOR = 'OpenBLAS'
    OMPI_MCA_plm = 'isolated'
    OMPI_MCA_btl = 'self'
  }
  post {
    always {
      juint 'CTestResults.xml'
    }

    failure {
      emailext(
	  subject: 'Failure in job ${currentBuild.fullDisplayName}',
	  recipientProviders: [[$class: 'CulpritsRecipientProvider']],
	  to: 'akantu-admins@akantu.ch',
	  attachLog: true,
          compressLog: true,
	  body: 'Something is wrong with ${env.BUILD_URL}')
    }
  }
}
