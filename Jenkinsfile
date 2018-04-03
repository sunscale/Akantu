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
        warnings(consoleParsers: [[parserName: 'GNU Make + GNU C Compiler (gcc)']])
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
        sh 'cd build/ && ctest -T test --no-compress-output || true'
	sh 'cp build/Testing/`head -n 1 < build/Testing/TAG`/Test.xml CTestResults.xml'
	// sh 'curl https://raw.githubusercontent.com/rpavlik/jenkins-ctest-plugin/master/ctest-to-junit.xsl -o ctest-to-junit.xsl'
	// sh 'xsltproc ctest-to-junit.xsl  build/Testing/`head -n 1 < build/Testing/TAG`/Test.xml > CTestResults.xml'
      }
    }
  }
  environment {
    BLA_VENDOR = 'OpenBLAS'
    OMPI_MCA_plm = 'isolated'
    OMPI_MCA_btl = 'tcp,self'
  }
  post {
    always {
      step([$class: 'XUnitBuilder',
         thresholds: [
             [$class: 'SkippedThreshold', failureThreshold: '0'],
             [$class: 'FailedThreshold', failureThreshold: '0']],
          tools: [[$class: 'CTestType', pattern: 'CTestResults.xml']]])
      step([$class: 'XUnitBuilder',
         thresholds: [
             [$class: 'SkippedThreshold', failureThreshold: '100'],
             [$class: 'FailedThreshold', failureThreshold: '0']],
          tools: [[$class: 'GoogleTestType', pattern: 'build/gtest_reports/**']]])
    }

    failure {
      emailext(
          body: '''${SCRIPT, template="groovy-html.template"}''',
	  mimeType: 'text/html',
          subject: "[Jenkins] ${currentBuild.fullDisplayName} Failed",
	  recipientProviders: [[$class: 'CulpritsRecipientProvider']],
	  to: 'akantu-admins@akantu.ch',
	  replyTo: 'akantu-admins@akantu.ch',
	  attachLog: true,
          compressLog: false)
    }
  }
}
