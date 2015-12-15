drop database if exists dream;
create database dream;
use dream;

create table drug_name
(
  id int primary key not null auto_increment,
  name varchar(10) not null
)
engine = InnoDB;

create table cell_name
(
  id int primary key not null auto_increment,
  name varchar(10) not null
)
engine = InnoDB;

create table drug_feature
(
  id int primary key not null auto_increment,
  feature varchar(10) not null
)
engine = InnoDB;

create table cell_feature
(
  id int primary key not null auto_increment,
  feature varchar(10) not null
)
engine = InnoDB;

create table drug_drug_feature
(
  id int primary key not null auto_increment,
  feature varchar(10) not null
)
engine = InnoDB;

create table drug_cell_feature
(
  id int primary key not null auto_increment,
  feature varchar(10) not null
) engine = InnoDB;

create table drug_drug_cell_feature
(
  id int primary key not null auto_increment,
  feature varchar(10) not null
)
engine = InnoDB;

create table drug
(
  id int primary key not null auto_increment,
  drug_id int not null,
  feature_id int not null,
  value varchar(10),
   
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
  value varchar(10),
   
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
  value varchar(10),
   
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
  value varchar(10),
   
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
  value varchar(10) not null,
   
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

-- ----------------------------------------
-- Insert temp data to base tables
-- ----------------------------------------
insert into drug_name (name) 
values ('drugA'),('drugB'),('drugC');
insert into cell_name (name)
values ('cellA'),('cellB'),('cellC');
insert into drug_feature (feature)
values ('df1'),('df2'),('df3');
insert into cell_feature (feature)
values ('cf1'),('cf2'),('cf3');
insert into drug_drug_feature (feature)
values ('ddf1'),('ddf2'),('ddf3');
insert into drug_cell_feature (feature)
values ('dcf1'),('dcf2'),('dcf3');
insert into drug_drug_cell_feature (feature)
values ('ddcf1'),('ddcf2'),('ddcf3');

-- ----------------------------------------
-- Insert temp data to value tables
-- ----------------------------------------
insert into drug (drug_id,feature_id,value)
values (1,1,'d1f150'),(1,2,'d1f210'),(2,1,'d2f122'),(1,3,'d1f39'),(3,2,'d3f211');
insert into cell (cell_id,feature_id,value)
values (1,1,'c1f150'),(1,2,'c1f210'),(2,1,'c2f122'),(1,3,'c1f39'),(3,2,'c3f211');
insert into drug_drug (drugA_id,drugB_id,feature_id,value)
values (1,1,2,'d11f150'),(1,3,2,'d13f210'),(2,2,1,'d22f122'),(2,1,3,'d21f39'),(3,1,2,'d31f211');
insert into drug_cell (drug_id,cell_id,feature_id,value)
values (1,2,1,'d1c2f150'),(1,1,2,'d1c2f210'),(2,3,1,'d2c1f122'),(1,3,2,'d1c3f39'),(3,2,1,'d3c2f211');
