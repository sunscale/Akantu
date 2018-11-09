pipeline {
  parameters {string(defaultValue: '', description: 'api-token', name: 'API_TOKEN')
              string(defaultValue: '', description: 'buildable phid', name: 'TARGET_PHID')
              string(defaultValue: '', description: 'Commit id', name: 'COMMIT_ID')
              string(defaultValue: '', description: 'Diff id', name: 'DIFF_ID')
    }

  options {
    disableConcurrentBuilds()
  }

  agent {
    dockerfile {
      additionalBuildArgs '--tag akantu-environment'
    }
  }
  stages {
    stage('Configure') {
      steps {
        sh """
        env
        mkdir -p build
        cd build
        cmake -DAKANTU_COHESIVE_ELEMENT:BOOL=TRUE -DAKANTU_IMPLICIT:BOOL=TRUE -DAKANTU_PARALLEL:BOOL=TRUE -DAKANTU_PYTHON_INTERFACE:BOOL=TRUE -DAKANTU_TESTS:BOOL=TRUE ..
        """
      }
    }
    stage('Compile') {
      steps {
	sh 'make -C build/src || true'
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
        sh """
        rm -rf build/gtest_reports
        cd build/
        #source ./akantu_environement.sh
        ctest -T test --no-compress-output || true
        TAG=`head -n 1 < build/Testing/TAG`
        if [ -e build/Testing/${TAG}/Test.xml ]; then
	   cp build/Testing/${TAG}/Test.xml CTestResults.xml
        fi
        """
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
      archiveArtifacts artifacts: 'build/Testing/**', fingerprint: true
      createartifact()
    }

    success {
      send_fail_pass('pass')
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
      send_fail_pass('fail')
    }
  }
}

def send_fail_pass(state) {
    sh """
set +x
curl https://c4science.ch/api/harbormaster.sendmessage \
-d api.token=${API_TOKEN} \
-d buildTargetPHID=${TARGET_PHID} \
-d type=${state}
"""
}

def createartifact() {
    sh """ set +x
curl https://c4science.ch/api/harbormaster.createartifact \
-d api.token=${API_TOKEN} \
-d buildTargetPHID=${TARGET_PHID} \
-d artifactKey="Jenkins URI" \
-d artifactType=uri \
-d artifactData[uri]=${BUILD_URL} \
-d artifactData[name]="View Jenkins result" \
-d artifactData[ui.external]=1
"""
}
