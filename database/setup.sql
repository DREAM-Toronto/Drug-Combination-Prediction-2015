drop database if exists dream2015;
create database dream2015;
use dream2015;

create table drug_name
(
  id int primary key not null auto_increment,
  name varchar(25) unique not null
)
engine = InnoDB;

create table cell_name
(
  id int primary key not null auto_increment,
  name varchar(25) unique not null
)
engine = InnoDB;

create table drug_feature
(
  id int primary key not null auto_increment,
  feature varchar(50) unique not null
)
engine = InnoDB;

create table cell_feature
(
  id int primary key not null auto_increment,
  feature varchar(50) unique not null
)
engine = InnoDB;

create table drug_drug_feature
(
  id int primary key not null auto_increment,
  feature varchar(50) unique not null
)
engine = InnoDB;

create table drug_cell_feature
(
  id int primary key not null auto_increment,
  feature varchar(50) unique not null
) engine = InnoDB;

create table drug_drug_cell_feature
(
  id int primary key not null auto_increment,
  feature varchar(50) unique not null
)
engine = InnoDB;

create table drug
(
  id int primary key not null auto_increment,
  drug_id int not null,
  feature_id int not null,
  value decimal(20,10) not null,
   
  foreign key (drug_id)
  references drug_name (id)
  on delete no action
  on update no action,

  foreign key (feature_id)
  references drug_feature (id)
  on delete no action
  on update no action
) engine = InnoDB;

create table cell 
(
  id int primary key not null auto_increment,
  cell_id int not null,
  feature_id int not null,
  value decimal(20,10) not null,
   
  foreign key (cell_id)
  references cell_name (id)
  on delete no action
  on update no action,

  foreign key (feature_id)
  references cell_feature (id)
  on delete no action
  on update no action
) engine = InnoDB;

create table drug_drug
(
  id int primary key not null auto_increment,
  drugA_id int not null,
  drugB_id int not null,
  feature_id int not null,
  value decimal(20,10) not null,
   
  foreign key (drugA_id)
  references drug_name (id)
  on delete no action
  on update no action,

  foreign key (drugB_id)
  references drug_name (id)
  on delete no action
  on update no action,

  foreign key (feature_id)
  references drug_drug_feature (id)
  on delete no action
  on update no action
) engine = InnoDB;

create table drug_cell
(
  id int primary key not null auto_increment,
  drug_id int not null,
  cell_id int not null,
  feature_id int not null,
  value decimal(20,10) not null,
   
  foreign key (drug_id)
  references drug_name (id)
  on delete no action
  on update no action,

  foreign key (cell_id)
  references cell_name (id)
  on delete no action
  on update no action,

  foreign key (feature_id)
  references drug_cell_feature (id)
  on delete no action
  on update no action
) engine = InnoDB;

create table drug_drug_cell
(
  id int primary key not null auto_increment,
  drugA_id int not null,
  drugB_id int not null,
  cell_id int not null,
  feature_id int not null,
  value decimal(20,10) not null,
   
  foreign key (drugA_id)
  references drug_name (id)
  on delete no action
  on update no action,

  foreign key (drugB_id)
  references drug_name (id)
  on delete no action
  on update no action,

  foreign key (cell_id)
  references cell_name (id)
  on delete no action
  on update no action,

  foreign key (feature_id)
  references drug_drug_cell_feature (id)
  on delete no action
  on update no action
) engine = InnoDB;