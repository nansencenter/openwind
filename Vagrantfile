# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|

  config.vm.box = "ubuntu/xenial64"
  #config.vm.box_url = "https://app.vagrantup.com/bento/boxes/ubuntu-18.10"

  config.vm.define "openwind", primary: true do |openwind|
    openwind.vm.network :private_network, ip: "192.168.33.10"
    openwind.vm.network "forwarded_port", guest: 8888, host: 8888
  end


  config.vm.provider "virtualbox" do |v|
    v.memory = 5000
    v.cpus = 2
  end

  # If true, then any SSH connections made will enable agent forwarding.
  config.ssh.forward_agent = true
  config.ssh.forward_x11 = true

  config.vm.provision "ansible_local" do |ansible|
    ansible.playbook = "provisioning/site.yml"
    ansible.galaxy_role_file = 'provisioning/galaxy_requirements.yml'
    ansible.galaxy_command = 'ansible-galaxy install --role-file=%{role_file} --roles-path=%{roles_path} --force -c'
  end

end
